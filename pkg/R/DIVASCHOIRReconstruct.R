# Categorical reconstruction with fixed subject directions.
# Fit every retained component together; components add in logit space.
choir_reconstruct <- function(observed, V, order, ranks, max_iter = 1000L, tol = 1e-8) {
  # 1. Recover the observed indicators, even when extraction used MCA.
  if (identical(attr(observed, "divaschoir_method"), "mca")) {
    p <- attr(observed, "divaschoir_p_hat")
    z <- sweep(sweep(observed, 1L, sqrt(p), "*"), 1L, p, "+")
    if (anyNA(z) || any(abs(z - round(z)) > 1e-7))
      stop("MCA input cannot be recovered as one-hot observations.")
    z <- round(z)
    attr(z, "divaschoir_method") <- "nuclear_norm"
  } else z <- MatCenterDIVASCHOIR(observed, method = "nuclear_norm")
  metadata <- attr(z, "divaschoir_variable_levels")
  groups <- lapply(metadata, function(variable) match(variable$rows, rownames(z)))
  if (anyNA(unlist(groups)) || anyNA(z) || any(!z %in% c(0, 1)) ||
      any(vapply(groups, function(rows) any(colSums(z[rows, , drop = FALSE]) != 1), logical(1))))
    stop("Categorical reconstruction requires complete one-hot responses and category metadata.")
  if (any(rowSums(z) == 0)) stop("Remove unobserved levels before unpenalized categorical reconstruction.")
  n <- ncol(z); d <- nrow(z)
  if (is.null(V)) V <- matrix(0, n, 0L)
  if (!is.matrix(V) || nrow(V) != n || any(!is.finite(V)) ||
      length(order) != length(ranks) || anyNA(ranks) || any(!is.finite(ranks)) || any(ranks < 1L) ||
      any(ranks != as.integer(ranks)) || sum(ranks) != ncol(V))
    stop("Categorical reconstruction basis, component order and ranks do not match.")
  if (!is.null(rownames(V)) && !is.null(colnames(z)) &&
      !identical(rownames(V), colnames(z))) stop("Categorical responses and subject scores are misaligned.")
  # Center fixed scores so eta has zero subject mean; retain each component span.
  score_means <- colMeans(V)
  V <- sweep(V, 2L, score_means, "-")
  r <- ncol(V)
  if (r > 0L && qr(V)$rank < r) stop("Centered categorical score directions must be linearly independent.")
  mu <- log(rowMeans(z))
  for (rows in groups) mu[rows] <- mu[rows] - mean(mu[rows])
  baseline_theta <- matrix(mu, d, n)
  baseline_nll <- -sum(z * choir_log_softmax(baseline_theta, groups))

  # 2. Convex multinomial likelihood with within-variable identifiability.
  B <- matrix(0, d, r)
  status <- "intercept_only"
  if (r > 0L) {
    mu_variable <- CVXR::Variable(c(d, 1L))
    B_variable <- CVXR::Variable(c(d, r))
    theta <- mu_variable %*% matrix(1, 1L, n) + B_variable %*% t(V)
    group_matrix <- matrix(0, length(groups), d)
    loss <- 0
    for (j in seq_along(groups)) {
      rows <- groups[[j]]
      group_matrix[j, rows] <- 1
      logits <- theta[rows, , drop = FALSE]
      loss <- loss + sum(CVXR::log_sum_exp(logits, axis = 2L, keepdims = TRUE)) -
        sum(z[rows, , drop = FALSE] * logits)
    }
    problem <- CVXR::Problem(CVXR::Minimize(loss), list(
      group_matrix %*% mu_variable == 0, group_matrix %*% B_variable == 0))
    CVXR::psolve(problem, solver = "CLARABEL", num_iter = max_iter,
                 abstol = tol, reltol = tol, feastol = tol)
    status <- CVXR::status(problem)
    if (!identical(status, "optimal"))
      stop("Categorical reconstruction did not converge; CVXR status: ", status)
    mu <- as.numeric(CVXR::value(mu_variable))
    B <- matrix(CVXR::value(B_variable), d, r)
    for (rows in groups) {
      mu[rows] <- mu[rows] - mean(mu[rows])
      B[rows, ] <- sweep(B[rows, , drop = FALSE], 2L, colMeans(B[rows, , drop = FALSE]), "-")
    }
  }

  # 3. Reconstruct additive natural-parameter components and probabilities.
  dimnames(B) <- list(rownames(z), colnames(V))
  eta <- B %*% t(V)
  dimnames(eta) <- dimnames(z)
  theta <- sweep(eta, 1L, mu, "+")
  log_probs <- choir_log_softmax(theta, groups)
  near_boundary <- max(abs(theta)) > 30 || any(exp(log_probs) < 1e-10)
  if (near_boundary) warning("Categorical reconstruction has extreme logits; a finite unpenalized MLE may not exist.")
  blocks <- loadings <- list()
  offset <- 0L
  for (j in seq_along(order)) {
    cols <- offset + seq_len(ranks[j])
    key <- as.character(order[j])
    loadings[[key]] <- B[, cols, drop = FALSE]
    blocks[[key]] <- B[, cols, drop = FALSE] %*% t(V[, cols, drop = FALSE])
    dimnames(blocks[[key]]) <- dimnames(z)
    offset <- offset + ranks[j]
  }
  list(matBlockMap = blocks, matLoadingMap = loadings, mu = setNames(mu, rownames(z)),
       eta = eta, theta = theta, probabilities = exp(log_probs), residual = z - exp(log_probs),
       basis = V, score_means = score_means, nll = -sum(z * log_probs),
       baseline_nll = baseline_nll, status = status, near_boundary = near_boundary,
       method = "multinomial_fixed_scores")
}
