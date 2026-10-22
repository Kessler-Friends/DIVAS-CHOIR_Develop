#' Signal Extraction for DIVASCHOIR Categorical Blocks
#'
#' Extracts categorical signal spaces from already one-hot encoded, centered,
#' and scaled categorical data blocks.
#'
#' @param datablock A list of centered and scaled categorical matrices with
#'   categories/features in rows and samples in columns.
#' @param nsim An integer specifying the number of bootstrap samples.
#' @param iplot A logical to indicate whether plots should be generated for
#'   visualizing singular values.
#' @param cull A numeric value for the culling parameter to adjust signal rank.
#' @param seed Optional. An integer to set the seed for the random number
#'   generator.
#' @param fixed_rank Optional. A positive integer, or one positive integer per
#'   categorical block, specifying the initial categorical signal rank before
#'   any optional null-angle culling step.
#' @param lambda Optional positive nuclear-norm penalty, one per block or a
#'   shared value. Singular values are soft-thresholded by this value.
#' @param cull_rank Apply rank culling using simulated categorical null angles.
#'   Defaults to FALSE, retaining all initial signal directions.
#'
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
#'       from 1,000 null MCA leading directions to one fixed rank-r subspace
#'       estimated from an independent null draw under theta = mu.}
#'   }
#'
#' @export
DIVASCHOIRSignalExtract <- function(
    datablock, nsim = 400,
    iplot = FALSE, cull = 0.382, seed = NULL, fixed_rank = NULL,
    lambda = NULL, cull_rank = FALSE
) {

  if (!is.null(seed)) {
    set.seed(seed)
  }

  fitted_category_probs <- function(signal_hat, p_hat, variable_levels) {
    lapply(variable_levels, function(variable) {
      rows <- variable$rows
      levels <- variable$levels
      row_idx <- match(rows, rownames(signal_hat))
      if (anyNA(row_idx)) {
        stop("Could not match categorical levels back to MCA rows for variable ", variable$variable, ".")
      }

      p_variable <- p_hat[row_idx]
      mu_clr <- log(p_variable) - mean(log(p_variable))
      # Map the fitted MCA signal back to natural coordinates once.
      covariance <- diag(p_variable) - tcrossprod(p_variable)
      eig <- eigen(covariance, symmetric = TRUE)
      cutoff <- max(eig$values) * 1e-10
      inverse_sqrt <- ifelse(eig$values > cutoff,
                             1 / sqrt(pmax(eig$values, cutoff)), 0)
      back_map <- (eig$vectors * rep(inverse_sqrt, each = length(p_variable))) %*%
        t(eig$vectors)
      eta <- back_map %*% signal_hat[row_idx, , drop = FALSE]
      eta <- sweep(eta, 2L, colMeans(eta), "-")
      vapply(seq_len(ncol(signal_hat)), function(sample_idx) {
        logits <- mu_clr + eta[, sample_idx]
        exp_logits <- exp(logits - max(logits))
        exp_logits / sum(exp_logits)
      }, numeric(length(levels)))
    })
  }

  simulate_mca_multinomial <- function(probabilities, variable_levels, sample_ids) {
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
    MatCenterDIVASCHOIR(sim_block, allow_empty_levels = TRUE)
  }

  # Check input dimensions
  if(!methods::is(datablock, "list")) stop("Input datablock has to be a list with length >= 1.")
  if(length(datablock) < 1) stop("Input datablock has to be a list with length >= 1.")
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
  if (!is.null(lambda)) {
    if (!is.null(fixed_rank)) stop("`lambda` and `fixed_rank` cannot be used together.")
    if (length(lambda) == 1L) lambda <- rep(lambda, nb)
    if (length(lambda) != nb || anyNA(lambda) ||
        any(!is.finite(lambda)) || any(lambda <= 0)) {
      stop("`lambda` must contain one positive finite value per categorical block.")
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

    # Initial categorical signal rank by truncated SVD.
    sv <- svd(datablockc)

    rHat <- if (is.null(lambda)) {
      min(if (is.null(fixed_rank)) 5L else fixed_rank[ib], length(sv$d))
    } else sum(sv$d > lambda[ib])
    if (rHat == 0L) {
      stop("Nuclear penalty removed every categorical direction in block ",
           dataname[ib], "; use a smaller `categorical_lambda_grid`.")
    }
    idx <- seq_len(rHat)

    U_hat <- sv$u[, idx, drop = FALSE]
    d_hat <- if (is.null(lambda)) sv$d[idx] else sv$d[idx] - lambda[ib]
    V_hat <- sv$v[, idx, drop = FALSE]

    signal_hat <- U_hat %*%
      diag(d_hat, nrow = rHat) %*%
      t(V_hat)
    rownames(signal_hat) <- rownames(datablockc)
    colnames(signal_hat) <- colnames(datablockc)

    EHat <- datablockc - signal_hat
    singValsTilde <- d_hat

    # Under theta = mu, every sample has the same category probabilities.
    # A zero MCA signal leaves only mu_clr in the softmax logits.
    null_signal <- matrix(0, d, n, dimnames = dimnames(datablockc))
    probabilities_null <- fitted_category_probs(null_signal, p_hat, variable_levels)
    randAngleCache <- numeric(1000L)
    randAngleCacheLoad <- numeric(1000L)
    # Estimate the reference subspaces once and hold them fixed, as in
    # randDirAngleMJ. Categorical null fits supply the random directions.
    null_reference <- svd(simulate_mca_multinomial(
      probabilities_null, variable_levels, colnames(datablockc)
    ), nu = rHat, nv = rHat)
    cat("Simulating 1,000 categorical null directions against fixed subspaces (theta = mu).\n")
    for (s in seq_along(randAngleCache)) {
      null_fit1 <- svd(simulate_mca_multinomial(
        probabilities_null, variable_levels, colnames(datablockc)
      ), nu = 1L, nv = 1L)
      # Angle from each leading unit direction to the fixed rank-r subspace:
      # acos of the orthogonal projection length.
      randAngleCache[s] <- acosd(sqrt(sum(
        crossprod(null_reference$v, null_fit1$v)^2
      )))
      randAngleCacheLoad[s] <- acosd(sqrt(sum(
        crossprod(null_reference$u, null_fit1$u)^2
      )))
    }
    nullAngleCaches[[ib]] <- list(subject = randAngleCache,
                                 feature = randAngleCacheLoad)
    randAngle <- stats::quantile(randAngleCache, 0.05)
    randAngleLoad <- stats::quantile(randAngleCacheLoad, 0.05)

    rSteps[[ib]] <- rHat

    # Categorical parametric bootstrap for perturbation-angle bounds.
    # Retain the initial rank unless culling is explicitly requested.
    PCAnglesCacheFullBoot <- matrix(90, nsim, rHat)
    PCAnglesCacheFullBootLoad <- matrix(90, nsim, rHat)

    cat("Progress Through Bootstrapped Matrices:\n")
    cat(paste0("\n", strrep(".", nsim), "\n\n"))

    # The fitted probabilities are fixed; only categorical draws vary.
    probabilities_hat <- fitted_category_probs(signal_hat, p_hat, variable_levels)
    for (s in seq_len(nsim)) {
      randX <- simulate_mca_multinomial(probabilities_hat, variable_levels,
                                       colnames(datablockc))
      randU <- U_hat
      randV <- V_hat
      svdRand <- svd(randX)
      randUHat <- svdRand$u[, idx, drop = FALSE]
      randVHat <- svdRand$v[, idx, drop = FALSE]

      for (j in seq_len(rHat)) {
        svd_randV <- svd(t(randV) %*% randVHat[, seq_len(j), drop = FALSE])
        PCAnglesCacheFullBoot[s, j] <- acosd(min(svd_randV$d))

        svd_randU <- svd(t(randU) %*% randUHat[, seq_len(j), drop = FALSE])
        PCAnglesCacheFullBootLoad[s, j] <- acosd(min(svd_randU$d))
      }

      cat("\b|\n")
    }

    randAngle <- as.numeric(randAngle)
    randAngleLoad <- as.numeric(randAngleLoad)
    cull <- as.numeric(cull)

    subject_q95 <- apply(PCAnglesCacheFullBoot, 2L, stats::quantile, probs = 0.95)
    feature_q95 <- apply(PCAnglesCacheFullBootLoad, 2L, stats::quantile, probs = 0.95)
    subject_cutoff <- randAngle * cull
    feature_cutoff <- randAngleLoad * cull
    validPC <- if (cull_rank) subject_q95 < subject_cutoff else rep(TRUE, rHat)
    if (cull_rank) validPC[1L] <- TRUE
    retained <- which(validPC)
    rBar <- length(retained)
    phiBar <- stats::quantile(PCAnglesCacheFullBoot[, rBar], 0.95)
    psiBar <- stats::quantile(PCAnglesCacheFullBootLoad[, rBar], 0.95)
    rankDiagnostics[[ib]] <- data.frame(
      direction = seq_len(rHat), subject_q95_deg = as.numeric(subject_q95),
      feature_q95_deg = as.numeric(feature_q95),
      subject_random_deg = randAngle, feature_random_deg = randAngleLoad,
      subject_cutoff_deg = subject_cutoff, feature_cutoff_deg = feature_cutoff,
      retained = validPC
    )

    cat(sprintf("Culled Rank is %d.\n", rBar))

    VVHatCacheBar <- vector("list", nsim)
    UUHatCacheBar <- vector("list", nsim)
    singValsTildeBar <- singValsTilde[retained]
    probabilities_bar <- probabilities_hat
    if (rBar < rHat) {
      signal_bar <- U_hat[, retained, drop = FALSE] %*%
        diag(singValsTildeBar, nrow = length(singValsTildeBar)) %*%
        t(V_hat[, retained, drop = FALSE])
      rownames(signal_bar) <- rownames(datablockc)
      colnames(signal_bar) <- colnames(datablockc)
      probabilities_bar <- fitted_category_probs(signal_bar, p_hat, variable_levels)
    }
    for (s in seq_len(nsim)) {
      randX <- simulate_mca_multinomial(probabilities_bar, variable_levels,
                                       colnames(datablockc))
      randU <- U_hat[, retained, drop = FALSE]
      randV <- V_hat[, retained, drop = FALSE]
      svdRand <- svd(randX)
      randUHat <- svdRand$u[, seq_len(rBar), drop = FALSE]
      randVHat <- svdRand$v[, seq_len(rBar), drop = FALSE]
      VVHatCacheBar[[s]] <- t(randV) %*% randVHat
      UUHatCacheBar[[s]] <- t(randU) %*% randUHat
    }

    VBars[[ib]] <- V_hat[, validPC, drop = FALSE]
    UBars[[ib]] <- U_hat[, validPC, drop = FALSE]
    phiBars[ib] <- phiBar
    psiBars[ib] <- psiBar
    rBars[ib] <- rBar
    EHats[[ib]] <- EHat
    singVals[[ib]] <- sv$d
    singValsHat[[ib]] <- c(d_hat, rep(0, length(sv$d) - rHat))
    VVHatCacheBars[[ib]] <- VVHatCacheBar
    UUHatCacheBars[[ib]] <- UUHatCacheBar
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

  # Return the results as a list
  return(list(
    VBars = VBars,
    UBars = UBars,
    phiBars = phiBars,
    psiBars = psiBars,
    EHats = EHats,
    rBars = rBars,
    singVals = singVals,
    singValsHat = singValsHat,
    rSteps = rSteps,
    VVHatCacheBars = VVHatCacheBars,
    UUHatCacheBars = UUHatCacheBars,
    rankDiagnostics = rankDiagnostics,
    nullAngleCaches = nullAngleCaches
  ))
}
