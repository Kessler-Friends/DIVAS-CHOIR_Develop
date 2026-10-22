# Method-specific CHOIR estimators. Preprocessing lives in SignalExtractUtils.R;
# null simulations and bootstrap orchestration live in DIVASCHOIRSignalExtract.R.

# MCA: truncate the SVD of the centered/scaled indicator matrix.
choir_estimate_mca <- function(x, rank) {
  decomposition <- svd(x)
  rank <- min(rank, length(decomposition$d))
  values <- decomposition$d[seq_len(rank)]
  idx <- seq_len(rank)
  u <- decomposition$u[, idx, drop = FALSE]
  v <- decomposition$v[, idx, drop = FALSE]
  signal <- (u * rep(values, each = nrow(x))) %*% t(v)
  dimnames(signal) <- dimnames(x)
  list(u = u, v = v, d = values, spectrum = decomposition$d,
       signal = signal, rank = rank, method = "mca", lambda = NULL)
}

# Nuclear: multinomial NLL + lambda * ||eta||_*, with an unpenalized intercept.
# rank is accepted for interface compatibility; lambda determines the fitted rank.
choir_estimate_nuclear_norm <- function(x, rank = NULL, lambda = NULL,
                                        max_iter = 1000L, tol = 1e-7,
                                        lambda_grid = NULL) {
  # Prepare one-hot observations and logit-centering constraints.
  z <- MatCenterDIVASCHOIR(x, method = "nuclear_norm")
  setup <- attr(z, "divaschoir_setup")
  groups <- setup$groups

  # Fit one penalty. A held-out response has zero counts in every category.
  fit_counts <- function(counts, penalty, initial_theta = NULL) {
    q <- nrow(counts); n <- ncol(counts)
    eta_variable <- CVXR::Variable(c(q, n))
    mu_variable <- CVXR::Variable(c(q, 1L))
    theta_expression <- eta_variable + mu_variable %*% matrix(1, 1L, n)
    group_matrix <- matrix(0, length(groups), q)
    loss <- 0
    # Add each variable's multinomial negative log likelihood: log-sum-exp
    # of its category logits minus the observed category's logit. Totals are
    # 1 for observed responses and 0 for held-out responses, which add no loss.
    # Also record category groups for the centering constraints below.
    for (j in seq_along(groups)) {
      rows <- groups[[j]]
      group_matrix[j, rows] <- 1
      logits <- theta_expression[rows, , drop = FALSE]
      totals <- matrix(colSums(counts[rows, , drop = FALSE]), 1L, n)
      loss <- loss + sum(totals * CVXR::log_sum_exp(logits, axis = 2L, keepdims = TRUE)) -
        sum(counts[rows, , drop = FALSE] * logits)
    }
    objective <- CVXR::Minimize(loss + penalty * CVXR::norm_nuc(eta_variable))
    constraints <- list(
      eta_variable %*% matrix(1, n, 1L) == 0, # Subject-centered latent effects.
      group_matrix %*% eta_variable == 0,     # Category-centered latent effects.
      group_matrix %*% mu_variable == 0       # Category-centered intercepts.
    )
    problem <- CVXR::Problem(objective, constraints)

    if (!is.null(initial_theta)) {
      CVXR::value(mu_variable) <- matrix(rowMeans(initial_theta), q, 1L)
      CVXR::value(eta_variable) <- choir_center_logits(initial_theta, groups)
    }
    CVXR::psolve(problem, solver = 'CLARABEL', warm_start = TRUE,
                 num_iter = max_iter, abstol = tol, reltol = tol, feastol = tol)
    converged <- identical(CVXR::status(problem), 'optimal')
    if (!converged) stop('Nuclear fit did not converge; CVXR status: ', CVXR::status(problem))
    eta <- choir_center_logits(CVXR::value(eta_variable), groups)
    mu <- as.numeric(CVXR::value(mu_variable))
    for (rows in groups) mu[rows] <- mu[rows] - mean(mu[rows])
    theta <- sweep(eta, 1L, mu, '+')
    dimnames(eta) <- dimnames(theta) <- dimnames(z)
    log_probs <- choir_log_softmax(theta, groups)
    for (key in c("divaschoir_variable_levels", "divaschoir_p_hat"))
      attr(eta, key) <- attr(z, key)
    decomposition <- svd(eta)
    keep <- which(decomposition$d > max(1, decomposition$d[1L]) * 1e-7)
    list(u = decomposition$u[, keep, drop = FALSE],
         v = decomposition$v[, keep, drop = FALSE], d = decomposition$d[keep],
         spectrum = decomposition$d, signal = eta, rank = length(keep),
         method = "nuclear_norm", lambda = penalty, mu = mu, theta = theta,
         probability_matrix = exp(log_probs), log_probs = log_probs,
         objective = -sum(counts * log_probs) + penalty * sum(decomposition$d),
         converged = converged)
  }

  # Five-fold response CV; reuse folds and warm-start along the penalty grid.
  select_lambda_cv <- function() {
    if (is.null(lambda_grid)) {
      baseline <- z
      for (rows in groups) baseline[rows, ] <- rowMeans(z[rows, , drop = FALSE])
      lambda_max <- max(svd(choir_center_logits(z - baseline, groups), nu = 0L, nv = 0L)$d)
      lambda_grid <- max(lambda_max, 1e-6) * c(1.05, 0.75, 0.45, 0.25, 0.1)
    }
    grid <- sort(unique(lambda_grid), decreasing = TRUE)
    selection <- data.frame(lambda = grid, cv_loss = NA_real_)
    # Preserve caller RNG state; use labels only for scoring, never assignment.
    had_seed <- exists('.Random.seed', envir = .GlobalEnv, inherits = FALSE)
    if (had_seed) saved_seed <- get('.Random.seed', envir = .GlobalEnv)
    on.exit(if (had_seed) assign('.Random.seed', saved_seed, envir = .GlobalEnv)
            else if (exists('.Random.seed', envir = .GlobalEnv, inherits = FALSE))
              rm('.Random.seed', envir = .GlobalEnv), add = TRUE)
    set.seed(20261019)
    J <- length(groups); n <- ncol(z)
    if (J < 2L || n < 2L || J*n < 5L)
      stop('Five-fold response CV needs at least two variables, two subjects, and five responses; supply lambda explicitly.')
    repeat {
      folds <- matrix(sample(rep(1:5, length.out = n*J)), n, J)
      if (all(vapply(1:5, function(k) all(rowSums(folds != k) > 0) &&
                     all(colSums(folds != k) > 0), logical(1)))) break
    }
    starts <- vector('list', 5L)
    fold_scores <- matrix(NA_real_, length(grid), 5L)
    for (k in seq_along(grid)) {
      if (isTRUE(getOption("divaschoir.cv_progress", FALSE))) {
        cat(sprintf("Five-fold CV candidate %d/%d: lambda = %.8g\n",
                    k, length(grid), grid[k]))
        flush(stdout())
      }
      score <- 0
      for (fold in 1:5) {
        training <- z
        held <- matrix(FALSE, nrow(z), n)
        for (j in seq_len(J)) held[groups[[j]], folds[,j] == fold] <- TRUE
        training[held] <- 0
        held_out <- fit_counts(training, grid[k], starts[[fold]])
        starts[[fold]] <- held_out$theta
        fold_scores[k,fold] <- -sum(z[held] * held_out$log_probs[held])
        score <- score + fold_scores[k,fold]
      }
      selection$cv_loss[k] <- score / (n*J)
    }
    best <- which.min(selection$cv_loss)
    fit <- fit_counts(z, grid[best])
    fit$selection <- selection
    fit$selection_method <- "five_fold_response_cv"
    fit$cv_folds <- folds
    fit$cv_fold_scores <- fold_scores
    fit$cv_seed <- 20261019L
    fit
  }

  # Bootstrap/null refits pass a fixed lambda and skip CV.
  if (is.null(lambda)) return(select_lambda_cv())
  fit_counts(z, lambda)
}
