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
#'
#' @return A list containing:
#'   \describe{
#'     \item{VBars}{List of sample-space right singular vector matrices.}
#'     \item{UBars}{List of category-space left singular vector matrices.}
#'     \item{phiBars}{Vector of adjusted perturbation angles.}
#'     \item{psiBars}{Vector of loadings perturbation angles.}
#'     \item{EHats}{List of residual matrices after rank-5 reconstruction.}
#'     \item{rBars}{Vector of retained ranks for each categorical block.}
#'     \item{singVals}{List of singular values for each categorical block.}
#'     \item{singValsHat}{List of rank-truncated singular values.}
#'     \item{rSteps}{List of retained rank values.}
#'     \item{VVHatCacheBars}{List of empty bootstrap cache placeholders.}
#'     \item{UUHatCacheBars}{List of empty bootstrap cache placeholders.}
#'   }
#'
#' @export
DIVASCHOIRSignalExtract <- function(
    datablock, nsim = 400,
    iplot = FALSE, cull = 0.382, seed = NULL
) {

  if (!is.null(seed)) {
    set.seed(seed)
  }

  softmax_stable <- function(x) {
    x <- x - max(x)
    exp_x <- exp(x)
    exp_x / sum(exp_x)
  }

  simulate_mca_multinomial <- function(theta_hat, p_hat, variable_levels) {
    n <- ncol(theta_hat)
    sim_block <- vector("list", length(variable_levels))
    names(sim_block) <- vapply(variable_levels, `[[`, character(1), "variable")

    for (variable in variable_levels) {
      rows <- variable$rows
      levels <- variable$levels
      row_idx <- match(rows, rownames(theta_hat))
      if (anyNA(row_idx)) {
        stop("Could not match categorical levels back to MCA rows for variable ", variable$variable, ".")
      }

      sampled <- character(n)
      for (sample_idx in seq_len(n)) {
        eta <- log(p_hat[row_idx]) + sqrt(p_hat[row_idx]) * theta_hat[row_idx, sample_idx]
        probs <- softmax_stable(eta)
        sampled[sample_idx] <- sample(levels, size = 1L, prob = probs)
      }
      sim_block[[variable$variable]] <- factor(sampled, levels = levels)
    }

    sim_block <- as.data.frame(sim_block, stringsAsFactors = TRUE)
    rownames(sim_block) <- colnames(theta_hat)
    MatCenterDIVASCHOIR(sim_block)
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

    # Fixed initial categorical signal rank by truncated SVD.
    sv <- svd(datablockc)

    rHat <- min(5L, length(sv$d))
    idx <- seq_len(rHat)

    U_hat <- sv$u[, idx, drop = FALSE]
    d_hat <- sv$d[idx]
    V_hat <- sv$v[, idx, drop = FALSE]

    signal_hat <- U_hat %*%
      diag(d_hat, nrow = rHat) %*%
      t(V_hat)
    rownames(signal_hat) <- rownames(datablockc)
    colnames(signal_hat) <- colnames(datablockc)

    EHat <- datablockc - signal_hat
    singValsTilde <- d_hat

    randAngleCache <- randDirAngleMJ(n, rHat, 1000)
    randAngleCacheLoad <- randDirAngleMJ(d, rHat, 1000)
    randAngle <- stats::quantile(randAngleCache, 0.05)
    randAngleLoad <- stats::quantile(randAngleCacheLoad, 0.05)

    rSteps[[ib]] <- rHat

    # Bootstrap estimation for perturbation-angle bounds. The categorical
    # signal rank is fixed at rHat; bootstrapping should not cull dimensions.
    PCAnglesCacheFullBoot <- matrix(90, nsim, rHat)
    PCAnglesCacheFullBootLoad <- matrix(90, nsim, rHat)

    cat("Progress Through Bootstrapped Matrices:\n")
    cat(paste0("\n", strrep(".", nsim), "\n\n"))

    for (s in seq_len(nsim)) {
      randX <- simulate_mca_multinomial(signal_hat, p_hat, variable_levels)
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

    validPC <- rep(TRUE, rHat)
    rBar <- rHat
    phiBar <- stats::quantile(PCAnglesCacheFullBoot[, rBar], 0.95)
    psiBar <- stats::quantile(PCAnglesCacheFullBootLoad[, rBar], 0.95)

    cat(sprintf("Culled Rank is %d.\n", rBar))

    VVHatCacheBar <- vector("list", nsim)
    UUHatCacheBar <- vector("list", nsim)
    singValsTildeBar <- singValsTilde[seq_len(rBar)]
    for (s in seq_len(nsim)) {
      theta_bar <- U_hat[, seq_len(rBar), drop = FALSE] %*%
        diag(singValsTildeBar, nrow = length(singValsTildeBar)) %*%
        t(V_hat[, seq_len(rBar), drop = FALSE])
      rownames(theta_bar) <- rownames(datablockc)
      colnames(theta_bar) <- colnames(datablockc)
      randX <- simulate_mca_multinomial(theta_bar, p_hat, variable_levels)
      randU <- U_hat[, seq_len(rBar), drop = FALSE]
      randV <- V_hat[, seq_len(rBar), drop = FALSE]
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
    UUHatCacheBars = UUHatCacheBars
  ))
}
