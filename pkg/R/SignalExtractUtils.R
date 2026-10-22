#' Arccosine in Degrees
#'
#' The `acosd` function computes the inverse cosine (arccosine) of a numeric value and converts the result from radians to degrees.
#'
#' @param x A numeric value or vector. The input should be within the range [-1, 1], as values outside this range will produce NaN values.
#'
#' @return A numeric value or vector representing the arccosine of the input, expressed in degrees.
#'
acosd <- function(x) {
  x <- pmin(1, pmax(-1, x))
  return(acos(x) * 180 / pi)
}



#' Center a Matrix by Rows, Columns, or Both
#'
#' This function centers a matrix by subtracting the row means, column means, or both, depending on the specified centering options. It is useful in data preprocessing to normalize data for further analysis.
#'
#' @param X A numeric matrix to be centered.
#' @param iColCent A logical value indicating whether to center by columns. If `TRUE`, the function subtracts the row means from each column.
#' @param iRowCent A logical value indicating whether to center by rows. If `TRUE`, the function subtracts the column means from each row.
#'
#' @return A centered matrix with the specified adjustments applied. If both `iColCent` and `iRowCent` are `TRUE`, the matrix will be centered by both rows and columns.
#' @export
#'
MatCenterJP <- function(X, iColCent = F, iRowCent = F) {

  d <- nrow(X)
  n <- ncol(X)
  outMat <- X
  if (iColCent) {
    outMat <- sweep(outMat, 1, rowMeans(outMat), "-")
  }
  if (iRowCent) {
    outMat <- sweep(outMat, 2, colMeans(outMat), "-")
  }
  return(outMat)
}

#' Prepare a Categorical Block for the Selected CHOIR Estimator
#'
#' MCA returns centered, scaled indicators. Nuclear likelihood returns unchanged
#' one-hot observations; its intercept and latent-signal constraints are applied
#' during estimation, rather than by transforming the observations.
#'
#' For a block with `n` samples and total one-hot dimension `J`, this constructs
#' the matrix
#'
#' \deqn{\tilde{M}_k^T = D_{\hat{p}_k}^{-1/2}
#' ((X_k^{obs})^T - \hat{p}_k 1_n^T)}
#'
#' where `X_obs` is the `n x J` stacked one-hot matrix,
#' `p_hat = n^{-1} t(X_obs) 1_n`, and `D_p_hat = diag(p_hat)`.
#'
#' @param block A data.frame or matrix containing categorical variables measured
#'   across samples. Rows are samples and columns are categorical variables.
#' @param allow_empty_levels If TRUE, retain an unobserved factor level as a
#'   zero row. This is used for bootstrap samples so all draws have the same
#'   category rows as the observed block.
#'
#' @param method Either mca or nuclear_norm; selects preparation only.
#' @return A numeric matrix with categories/features in rows and samples in
#'   columns.
#'
#' @export
MatCenterDIVASCHOIR <- function(block, allow_empty_levels = FALSE,
                                method = c("mca", "nuclear_norm")) {
  method <- match.arg(method)
  prepare <- switch(method, mca = choir_prepare_mca, nuclear_norm = choir_prepare_nuclear)
  if (!is.null(attr(block, "divaschoir_variable_levels"))) {
    previous <- attr(block, "divaschoir_method")
    if (is.null(previous)) previous <- "mca"
    if (previous != method) stop("Prepared CHOIR block uses ", previous,
                                 "; pass raw data to change methods.")
    return(block)
  }
  if (!is.logical(allow_empty_levels) || length(allow_empty_levels) != 1L ||
      is.na(allow_empty_levels))
    stop("`allow_empty_levels` must be TRUE or FALSE.")
  block <- as.data.frame(block)
  block[] <- lapply(block, function(x) if (is.factor(x)) x else factor(x))
  indicator_list <- Map(function(x, variable_name) {
    mm <- stats::model.matrix(~ x - 1)
    level_names <- sub("^x", "", colnames(mm))
    colnames(mm) <- paste0(variable_name, level_names)
    mm
  }, block, names(block))
  indicator <- do.call(cbind, indicator_list)
  p_hat <- colMeans(indicator)
  if (any(p_hat == 0) && !allow_empty_levels)
    stop("A categorical level is unobserved; remove unused levels from the observed block.")
  x <- t(indicator)
  colnames(x) <- rownames(block)
  attr(x, "divaschoir_variable_levels") <- Map(function(variable_name, indicator_matrix) {
    list(variable = variable_name, levels = levels(block[[variable_name]]),
         rows = colnames(indicator_matrix))
  }, names(block), indicator_list)
  attr(x, "divaschoir_p_hat") <- p_hat
  prepare(x)
}

# MCA: center and scale the shared one-hot matrix.
choir_prepare_mca <- function(x) {
  p <- attr(x, "divaschoir_p_hat")
  x <- sweep(x, 1L, p, "-")
  x <- sweep(x, 1L, sqrt(ifelse(p > 0, p, 1)), "/")
  attr(x, "divaschoir_method") <- "mca"
  x
}

# Log-softmax across each variable's categories, separately for each subject.
choir_log_softmax <- function(theta, groups) {
  log_probs <- theta
  for (rows in groups) {
    logits <- theta[rows, , drop = FALSE]
    logits <- sweep(logits, 2L, apply(logits, 2L, max), "-")
    log_probs[rows, ] <- sweep(logits, 2L, log(colSums(exp(logits))), "-")
  }
  log_probs
}

# Center logits within each variable, optionally across subjects as well.
choir_center_logits <- function(x, groups, subjects = TRUE) {
  if (subjects) x <- sweep(x, 1L, rowMeans(x), "-")
  for (rows in groups) x[rows, ] <- sweep(x[rows, , drop = FALSE], 2L,
    colMeans(x[rows, , drop = FALSE]), "-")
  x
}

# Keep observations unchanged; initialize the constrained latent parameters.
choir_prepare_nuclear <- function(z) {
  groups <- lapply(attr(z, "divaschoir_variable_levels"),
                   function(variable) match(variable$rows, rownames(z)))
  mu <- log(pmax(rowMeans(z), 1e-8))
  for (rows in groups) mu[rows] <- mu[rows] - mean(mu[rows])
  attr(z, "divaschoir_setup") <- list(
    groups = groups,
    mu = mu,
    eta = matrix(0, nrow(z), ncol(z), dimnames = dimnames(z))
  )
  attr(z, "divaschoir_method") <- "nuclear_norm"
  z
}
