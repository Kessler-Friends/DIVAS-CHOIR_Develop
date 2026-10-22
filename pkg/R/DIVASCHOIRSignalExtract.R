#' Signal Extraction for DIVASCHOIR Categorical Blocks
#'
#' DIVASmain prepares observed categorical blocks before calling this function.
#' Direct calls may also supply raw data; new bootstrap draws are prepared here.
#' Extracts categorical signal spaces using method-specific preprocessing:
#' MCA centers/scales indicators; nuclear likelihood centers the latent signal.
#'
#' @param datablock A list of raw categorical data frames (subjects by variables),
#'   or tagged matrices prepared with MatCenterDIVASCHOIR for the selected method.
#' @param nsim An integer specifying the number of bootstrap samples.
#' @param iplot A logical to indicate whether plots should be generated for
#'   visualizing singular values.
#' @param cull A numeric value for the culling parameter to adjust signal rank; defaults to 0.5.
#' @param seed Optional. An integer to set the seed for the random number
#'   generator.
#' @param fixed_rank Optional. A positive integer, or one positive integer per
#'   categorical block. For nuclear norm, truncate the fitted latent signal to this rank.
#' @param cull_rank Apply rank culling using simulated categorical null angles.
#'   Defaults to FALSE, retaining all initial signal directions.
#'   No direction is forced to remain: if none pass, return rank zero and
#'   empty bases/caches, with 90-degree angle placeholders.
#' @param method CHOIR estimator: "mca" or "nuclear_norm". The same
#'   estimator is used for original, null and bootstrap data. Nuclear norm uses
#'   five-fold response CV to select lambda; refits retain the requested leading directions.
#' @param lambda_grid Optional nonnegative penalties for nuclear five-fold response CV selection.
#'   NULL uses factors 1.05, 0.75, 0.45, 0.25, and 0.1 times the spectral norm
#'   of the centered intercept-only score; ignored for MCA.
#'
#' @param simulate_null Compute null-angle diagnostics. Must be TRUE when cull_rank is TRUE.
#' @param bootstrap_tol Nuclear refit convergence tolerance. Original/CV fits use 1e-7.
#' @param null_nsim Positive integer number of categorical null simulations.
#'   Defaults to 1000; separate from the number of bootstrap samples.
#' @return A list containing:
#'   \describe{
#'     \item{VBars}{List of sample-space right singular vector matrices.}
#'     \item{UBars}{List of category-space left singular vector matrices.}
#'     \item{phiBars}{Vector of adjusted perturbation angles.}
#'     \item{psiBars}{Vector of loadings perturbation angles.}
#'     \item{EHats}{List of residual matrices after initial-rank reconstruction.}
#'     \item{rBars}{Vector of retained ranks for each categorical block.}
#'     \item{singVals}{List of singular values for each categorical block.}
#'     \item{singValsHat}{List of rank-truncated singular values.}
#'     \item{rSteps}{List of retained rank values.}
#'     \item{VVHatCacheBars}{List of empty bootstrap cache placeholders.}
#'     \item{UUHatCacheBars}{List of empty bootstrap cache placeholders.}
#'     \item{rankDiagnostics}{Per-block bootstrap angle quantiles and culling cutoffs.}
#'     \item{nullAngleCaches}{Per-block sample-space and category-space angles
#'       from null_nsim null leading directions to one fixed rank-r subspace
#'       estimated from an independent null draw under theta = mu.}
#'   }
#'
#' @export
DIVASCHOIRSignalExtract <- function(
    datablock, nsim = 400,
    iplot = FALSE, cull = 0.5, seed = NULL, fixed_rank = NULL,
    cull_rank = FALSE, method = c("mca", "nuclear_norm"), lambda_grid = NULL,
    simulate_null = TRUE, bootstrap_tol = 1e-7, null_nsim = 1000L
) {
  method <- match.arg(method)
  if (length(bootstrap_tol) != 1L || !is.finite(bootstrap_tol) || bootstrap_tol <= 0)
    stop("bootstrap_tol must be a finite positive number")
  if (length(null_nsim) != 1L || !is.numeric(null_nsim) || !is.finite(null_nsim) ||
      null_nsim < 1 || null_nsim != floor(null_nsim) || null_nsim > .Machine$integer.max)
    stop("null_nsim must be a positive integer")
  null_nsim <- as.integer(null_nsim)
  if (cull_rank && !isTRUE(simulate_null)) stop("Rank culling requires null simulations")
  prepare <- function(block, allow_empty_levels = FALSE)
    MatCenterDIVASCHOIR(block, allow_empty_levels, method)
  block_lambda <- NULL
  estimate <- if (method == "mca") choir_estimate_mca else function(x, rank = NULL) {
    fit <- choir_estimate_nuclear_norm(x, lambda = block_lambda, lambda_grid = lambda_grid,
      tol = if (is.null(block_lambda)) 1e-7 else bootstrap_tol)
    # Select lambda on the penalized likelihood, then truncate latent eta.
    # With no explicit initial rank, keep the penalty-selected initial rank.
    target <- if (is.null(block_lambda) && is.null(fixed_rank)) NULL else rank
    fit$untruncated_rank <- fit$rank
    if (!is.null(target) && fit$rank > target) {
      idx <- seq_len(target)
      fit$u <- fit$u[, idx, drop = FALSE]
      fit$v <- fit$v[, idx, drop = FALSE]
      fit$d <- fit$d[idx]
      fit$signal[] <- (fit$u * rep(fit$d, each = nrow(fit$u))) %*% t(fit$v)
      fit$rank <- length(idx)
      fit$theta <- sweep(fit$signal, 1L, fit$mu, "+")
      groups <- attr(x, "divaschoir_setup")$groups
      fit$log_probs <- choir_log_softmax(fit$theta, groups)
      fit$probability_matrix <- exp(fit$log_probs)
    }
    fit
  }

  if (!is.null(seed)) set.seed(seed)

  fitted_category_probs <- function(signal_hat, p_hat, variable_levels, mu = NULL) {
    lapply(variable_levels, function(variable) {
      row_idx <- match(variable$rows, rownames(signal_hat))
      if (anyNA(row_idx)) {
        stop("Could not match categorical levels back to MCA rows for variable ", variable$variable, ".")
      }

      eta <- signal_hat[row_idx, , drop = FALSE]
      baseline <- if (is.null(mu)) NULL else mu[row_idx]
      if (method == "mca") {
        p <- p_hat[row_idx]
        baseline <- log(p) - mean(log(p))
        # Map the fitted MCA signal back to natural coordinates.
        eig <- eigen(diag(p) - tcrossprod(p), symmetric = TRUE)
        cutoff <- max(eig$values) * 1e-10
        inverse_sqrt <- ifelse(eig$values > cutoff,
                               1 / sqrt(pmax(eig$values, cutoff)), 0)
        back_map <- (eig$vectors * rep(inverse_sqrt, each = length(p))) %*% t(eig$vectors)
        eta <- back_map %*% eta
        eta <- sweep(eta, 2L, colMeans(eta), "-")
      }
      if (is.null(baseline)) stop("Nuclear probabilities require the fitted intercept.")
      logits <- sweep(eta, 1L, baseline, "+")
      exp(choir_log_softmax(logits, list(seq_len(nrow(logits)))))
    })
  }

  simulate_categorical <- function(probabilities, variable_levels, sample_ids) {
    n <- ncol(probabilities[[1L]])
    sim_block <- vector("list", length(variable_levels))
    names(sim_block) <- vapply(variable_levels, `[[`, character(1), "variable")

    for (j in seq_along(variable_levels)) {
      variable <- variable_levels[[j]]
      levels <- variable$levels
      sampled <- character(n)
      for (sample_idx in seq_len(n)) {
        sampled[sample_idx] <- sample(levels, size = 1L,
                                      prob = probabilities[[j]][, sample_idx])
      }
      sim_block[[variable$variable]] <- factor(sampled, levels = levels)
    }

    sim_block <- as.data.frame(sim_block, stringsAsFactors = TRUE)
    rownames(sim_block) <- sample_ids
    prepare(sim_block, allow_empty_levels = TRUE)
  }

  if (!is.list(datablock) || !length(datablock)) stop("Supply a nonempty list of CHOIR blocks.")
  # DIVASmain supplies prepared observed blocks. Keep raw-input support for
  # direct calls, and prepare every new categorical bootstrap draw above.
  datablock <- lapply(datablock, function(block) {
    if (is.null(attr(block, "divaschoir_variable_levels"))) return(prepare(block))
    previous <- attr(block, "divaschoir_method")
    if (is.null(previous)) previous <- "mca"
    if (previous != method)
      stop("Prepared CHOIR block uses ", previous, "; pass raw data to change methods.")
    block
  })

  # Check input dimensions
  if(max(sapply(datablock, ncol)) != min(sapply(datablock, ncol))) stop("All data blocks have to have the same number of columns (samples).")

  # Check data block names
  nb <- length(datablock)
  dataname <- names(datablock)
  if(is.null(dataname)){
    warning("Input datablock is unnamed, generic names for categorical data blocks generated.")
    dataname <- paste0("CategoricalDatablock", 1:nb)
  }

  if (!is.null(fixed_rank)) {
    fixed_rank <- as.integer(fixed_rank)
    if (length(fixed_rank) == 1L) fixed_rank <- rep(fixed_rank, nb)
    if (length(fixed_rank) != nb || anyNA(fixed_rank) || any(fixed_rank < 1L)) {
      stop("`fixed_rank` must be a positive integer or one positive integer per categorical block.")
    }
  }
  # Initialize the output lists
  VBars <- vector("list", nb) # adjusted signal row spaces
  UBars <- vector("list", nb) # adjusted signal column spaces
  EHats <- vector("list", nb) # estimated residual matrices
  phiBars <- rep(90, nb) # adjusted perturbation angles
  psiBars <- rep(90, nb) # loadings perturbation angles
  rBars <- numeric(nb) # adjusted signal ranks
  singVals <- vector("list", nb) # singular values before truncation
  singValsHat <- vector("list", nb) # singular values after rank truncation
  rSteps <- vector("list", nb) # signal rank adjustment steps
  VVHatCacheBars <- vector("list", nb) # cached VVHat matrices for bootstrap samples
  UUHatCacheBars <- vector("list", nb) # cached UUHat matrices for bootstrap samples
  rankDiagnostics <- vector("list", nb)
  nullAngleCaches <- vector("list", nb)
  estimation <- vector("list", nb)
  analysisBlocks <- vector("list", nb)

  # Loop through each block
  for (ib in seq_len(nb)) {
    cat(sprintf("Categorical signal estimation for %s\n", dataname[ib]))

    datablockc <- datablock[[ib]]
    d <- nrow(datablockc)
    n <- ncol(datablockc)
    p_hat <- attr(datablockc, "divaschoir_p_hat")
    variable_levels <- attr(datablockc, "divaschoir_variable_levels")
    if (is.null(p_hat) || is.null(variable_levels)) {
      stop(
        "Categorical block ", dataname[ib],
        " is missing DIVASCHOIR level metadata needed for multinomial bootstrap. ",
        "Pass raw categorical blocks through DIVASmain/MatCenterDIVASCHOIR."
      )
    }

    cat("Categorical datablock dimensions:", d, "features;", n, "samples \n")

    block_lambda <- NULL
    initial_fit <- estimate(datablockc,
      rank = if (is.null(fixed_rank)) 5L else fixed_rank[ib])
    rHat <- initial_fit$rank
    block_lambda <- initial_fit$lambda

    U_hat <- initial_fit$u
    d_hat <- initial_fit$d
    V_hat <- initial_fit$v
    signal_hat <- initial_fit$signal

    EHat <- if (method == "mca") datablockc - signal_hat else datablockc - initial_fit$probability_matrix
    analysisBlocks[[ib]] <- if (method == "mca") datablockc else signal_hat

    estimation[[ib]] <- list(method = method, rank = rHat, untruncated_rank = initial_fit$untruncated_rank,
      rank_method = if (method == "nuclear_norm") "post_fit_svd_truncation" else "svd", lambda = initial_fit$lambda,
      mu = initial_fit$mu, converged = initial_fit$converged, objective = initial_fit$objective,
      selection = initial_fit$selection, selection_method = initial_fit$selection_method,
      cv_folds = initial_fit$cv_folds, cv_seed = initial_fit$cv_seed,
      null_nsim = if (isTRUE(simulate_null)) null_nsim else 0L)
    # A CV-selected intercept-only model has no subspaces to bootstrap.
    if (rHat == 0L) {
      UBars[[ib]] <- U_hat
      VBars[[ib]] <- V_hat
      EHats[[ib]] <- EHat
      singVals[[ib]] <- initial_fit$spectrum
      singValsHat[[ib]] <- rep(0, length(initial_fit$spectrum))
      rSteps[[ib]] <- 0L
      VVHatCacheBars[[ib]] <- UUHatCacheBars[[ib]] <- list()
      rankDiagnostics[[ib]] <- data.frame(direction = integer(), retained = logical())
      nullAngleCaches[[ib]] <- list(subject = numeric(), feature = numeric())
      estimation[[ib]]$angle_ranks <- estimation[[ib]]$retained_ranks <- integer()
      estimation[[ib]]$angle_lambdas <- estimation[[ib]]$retained_lambdas <- numeric()
      next
    }

    randAngle <- randAngleLoad <- NA_real_
    randAngleCache <- randAngleCacheLoad <- numeric()
    null_ranks <- integer()
    null_reference <- list(rank = NA_integer_)
    if (isTRUE(simulate_null)) {
      # Under theta = mu, every sample has the same category probabilities.
      # A zero MCA signal leaves only mu_clr in the softmax logits.
      null_signal <- matrix(0, d, n, dimnames = dimnames(datablockc))
      probabilities_null <- fitted_category_probs(null_signal, p_hat, variable_levels, initial_fit$mu)
      randAngleCache <- randAngleCacheLoad <- rep(NA_real_, null_nsim)
      null_ranks <- integer(null_nsim)
      # Estimate the reference subspaces once and hold them fixed, as in
      # randDirAngleMJ. Categorical null fits supply the random directions.
      null_reference <- estimate(simulate_categorical(
        probabilities_null, variable_levels, colnames(datablockc)
      ), rank = rHat)
      cat(sprintf("Simulating %d categorical null directions against fixed subspaces (theta = mu).\n", null_nsim))
      null_started <- proc.time()[[3L]]
      for (s in seq_along(randAngleCache)) {
        null_fit1 <- estimate(simulate_categorical(
          probabilities_null, variable_levels, colnames(datablockc)
        ), rank = 1L)
        null_ranks[s] <- null_fit1$rank
        if (s == 1L || s %% 10L == 0L || s == length(randAngleCache)) {
          cat(sprintf("Null simulation: %d/%d (%.1f seconds elapsed)\n",
                      s, length(randAngleCache), proc.time()[[3L]] - null_started))
          flush(stdout())
        }
        # Zero-rank null fits have no random direction; leave their angle undefined.
        if (null_reference$rank > 0L && null_fit1$rank > 0L) {
          reference_idx <- seq_len(min(rHat, null_reference$rank))
          randAngleCache[s] <- acosd(sqrt(sum(crossprod(
            null_reference$v[, reference_idx, drop = FALSE], null_fit1$v[, 1L])^2)))
          randAngleCacheLoad[s] <- acosd(sqrt(sum(crossprod(
            null_reference$u[, reference_idx, drop = FALSE], null_fit1$u[, 1L])^2)))
        }
      }
      nullAngleCaches[[ib]] <- list(subject = randAngleCache,
                                   feature = randAngleCacheLoad)
      randAngle <- stats::quantile(randAngleCache, 0.05, na.rm = TRUE)
      randAngleLoad <- stats::quantile(randAngleCacheLoad, 0.05, na.rm = TRUE)

      if (cull_rank && (!is.finite(randAngle) || !is.finite(randAngleLoad)))
        stop("Null fits have no signal directions at the fixed nuclear penalty; null-angle culling is unavailable.")
    } else {
      cat("Skipping categorical null simulations (rank culling disabled).\n")
      nullAngleCaches[[ib]] <- list(subject = numeric(), feature = numeric())
    }
    rSteps[[ib]] <- rHat

    # Categorical parametric bootstrap for perturbation-angle bounds.
    # Retain the initial rank unless culling is explicitly requested.
    PCAnglesCacheFullBoot <- matrix(90, nsim, rHat)
    PCAnglesCacheFullBootLoad <- matrix(90, nsim, rHat)


    # The fitted probabilities are fixed; only categorical draws vary.
    probabilities_hat <- fitted_category_probs(signal_hat, p_hat, variable_levels, initial_fit$mu)
    angle_lambdas <- rep(NA_real_, nsim)
    angle_ranks <- integer(nsim)
    for (s in seq_len(nsim)) {
      randX <- simulate_categorical(probabilities_hat, variable_levels,
                                       colnames(datablockc))
      svdRand <- estimate(randX, rank = rHat)
      if (s %% 10L == 0L) { cat(sprintf("Angle bootstrap: %d/%d\n", s, nsim)); flush(stdout()) }
      angle_ranks[s] <- svdRand$rank
      if (!is.null(svdRand$lambda)) angle_lambdas[s] <- svdRand$lambda
      available <- seq_len(min(rHat, svdRand$rank))
      randUHat <- svdRand$u[, available, drop = FALSE]
      randVHat <- svdRand$v[, available, drop = FALSE]

      # Missing bootstrap directions keep their initialized 90-degree bound.
      for (j in available) {
        svd_randV <- svd(crossprod(V_hat, randVHat[, seq_len(j), drop = FALSE]))
        PCAnglesCacheFullBoot[s, j] <- acosd(min(svd_randV$d))

        svd_randU <- svd(crossprod(U_hat, randUHat[, seq_len(j), drop = FALSE]))
        PCAnglesCacheFullBootLoad[s, j] <- acosd(min(svd_randU$d))
      }

    }

    randAngle <- as.numeric(randAngle)
    randAngleLoad <- as.numeric(randAngleLoad)
    cull <- as.numeric(cull)

    subject_q95 <- apply(PCAnglesCacheFullBoot, 2L, stats::quantile, probs = 0.95)
    feature_q95 <- apply(PCAnglesCacheFullBootLoad, 2L, stats::quantile, probs = 0.95)
    subject_cutoff <- randAngle * cull
    feature_cutoff <- randAngleLoad * cull
    validPC <- if (cull_rank) subject_q95 < subject_cutoff else rep(TRUE, rHat)
    retained <- which(validPC)
    rBar <- length(retained)
    # Rank zero has no angle column to summarize. Use the same empty-space
    # convention as an intercept-only initial fit.
    phiBar <- if (rBar > 0L) stats::quantile(PCAnglesCacheFullBoot[, rBar], 0.95) else 90
    psiBar <- if (rBar > 0L) stats::quantile(PCAnglesCacheFullBootLoad[, rBar], 0.95) else 90
    rankDiagnostics[[ib]] <- data.frame(
      direction = seq_len(rHat), subject_q95_deg = as.numeric(subject_q95),
      feature_q95_deg = as.numeric(feature_q95),
      subject_random_deg = randAngle, feature_random_deg = randAngleLoad,
      subject_cutoff_deg = subject_cutoff, feature_cutoff_deg = feature_cutoff,
      retained = validPC
    )

    cat(sprintf("Culled Rank is %d.\n", rBar))

    # There are no retained subspaces to bootstrap if all directions fail.
    retained_nsim <- if (rBar > 0L) nsim else 0L
    VVHatCacheBar <- vector("list", retained_nsim)
    UUHatCacheBar <- vector("list", retained_nsim)
    singValsTildeBar <- d_hat[retained]
    probabilities_bar <- probabilities_hat
    if (rBar > 0L && rBar < rHat) {
      signal_bar <- U_hat[, retained, drop = FALSE] %*%
        diag(singValsTildeBar, nrow = length(singValsTildeBar)) %*%
        t(V_hat[, retained, drop = FALSE])
      rownames(signal_bar) <- rownames(datablockc)
      colnames(signal_bar) <- colnames(datablockc)
      probabilities_bar <- fitted_category_probs(signal_bar, p_hat, variable_levels, initial_fit$mu)
    }
    retained_lambdas <- rep(NA_real_, retained_nsim)
    retained_ranks <- integer(retained_nsim)
    for (s in seq_len(retained_nsim)) {
      randX <- simulate_categorical(probabilities_bar, variable_levels,
                                       colnames(datablockc))
      svdRand <- estimate(randX, rank = rBar)
      if (s %% 10L == 0L) { cat(sprintf("Retained bootstrap: %d/%d\n", s, nsim)); flush(stdout()) }
      retained_ranks[s] <- svdRand$rank
      if (!is.null(svdRand$lambda)) retained_lambdas[s] <- svdRand$lambda
      available <- seq_len(min(rBar, svdRand$rank))
      # Cache zeros for missing directions rather than arbitrary zero-singular-value bases.
      VVHatCacheBar[[s]] <- UUHatCacheBar[[s]] <- matrix(0, rBar, rBar)
      if (length(available)) {
        VVHatCacheBar[[s]][, available] <- crossprod(V_hat[, retained, drop = FALSE],
                                                   svdRand$v[, available, drop = FALSE])
        UUHatCacheBar[[s]][, available] <- crossprod(U_hat[, retained, drop = FALSE],
                                                   svdRand$u[, available, drop = FALSE])
      }
    }

    VBars[[ib]] <- V_hat[, validPC, drop = FALSE]
    UBars[[ib]] <- U_hat[, validPC, drop = FALSE]
    phiBars[ib] <- phiBar
    psiBars[ib] <- psiBar
    rBars[ib] <- rBar
    EHats[[ib]] <- EHat
    singVals[[ib]] <- initial_fit$spectrum
    singValsHat[[ib]] <- c(d_hat, rep(0, length(initial_fit$spectrum) - rHat))
    VVHatCacheBars[[ib]] <- VVHatCacheBar
    UUHatCacheBars[[ib]] <- UUHatCacheBar
    estimation[[ib]]$bootstrap_tol <- bootstrap_tol
    estimation[[ib]]$simulate_null <- simulate_null
    estimation[[ib]]$angle_ranks <- angle_ranks
    estimation[[ib]]$angle_lambdas <- angle_lambdas
    estimation[[ib]]$retained_ranks <- retained_ranks
    estimation[[ib]]$retained_lambdas <- retained_lambdas
    estimation[[ib]]$null_ranks <- null_ranks
    estimation[[ib]]$null_reference_rank <- null_reference$rank
    estimation[[ib]]$null_valid_draws <- sum(is.finite(randAngleCache))
  }

  # Plotting (if iplot is TRUE)
  if (iplot) {
    for (ib in seq_len(nb)) {
      singValsI <- singVals[[ib]]
      singValsHatI <- singValsHat[[ib]]
      rBar <- rBars[ib]
      matName <- dataname[ib]

      plot(seq_along(singValsI), singValsI, type = "b", col = "blue", pch = 16,
           xlab = "Index", ylab = "Singular Value",
           main = paste0(matName, " Categorical Singular Values"))
      graphics::points(seq_along(singValsHatI), singValsHatI, type = "b", col = "red", pch = 4)
      graphics::abline(v = rBar, col = "black", lty = 3)
    }
  }

  list(
    VBars = VBars, UBars = UBars, phiBars = phiBars, psiBars = psiBars,
    EHats = EHats, rBars = rBars, singVals = singVals, singValsHat = singValsHat,
    rSteps = rSteps, VVHatCacheBars = VVHatCacheBars, UUHatCacheBars = UUHatCacheBars,
    rankDiagnostics = rankDiagnostics, nullAngleCaches = nullAngleCaches,
    estimation = estimation, preparedBlocks = datablock, analysisBlocks = analysisBlocks
  )
}
