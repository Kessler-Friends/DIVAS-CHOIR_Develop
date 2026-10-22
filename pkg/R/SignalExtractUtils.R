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



#' Construct MCA-Style SVD Input for Categorical DIVAS Blocks
#'
#' `MatCenterDIVASCHOIR` converts a categorical data block into the centered
#' and scaled indicator matrix used for MCA-style singular value decomposition.
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
#' @return A numeric matrix with categories/features in rows and samples in
#'   columns.
#'
#' @export
MatCenterDIVASCHOIR <- function(block, allow_empty_levels = FALSE) {
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

  svd_input <- sweep(t(indicator), 1, p_hat, "-")
  # An absent bootstrap category is already a zero row; divide it by one.
  svd_input <- sweep(svd_input, 1, sqrt(ifelse(p_hat > 0, p_hat, 1)), "/")

  colnames(svd_input) <- rownames(block)
  attr(svd_input, "divaschoir_variable_levels") <- Map(function(variable_name, indicator_matrix) {
    list(
      variable = variable_name,
      levels = levels(block[[variable_name]]),
      rows = colnames(indicator_matrix)
    )
  }, names(block), indicator_list)
  attr(svd_input, "divaschoir_p_hat") <- p_hat
  svd_input
}
