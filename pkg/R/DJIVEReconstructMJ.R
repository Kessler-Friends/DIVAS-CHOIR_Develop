#' MatReconstructMJ - Reconstruct joint block matrices and their loadings from data
#'
#' @param X d x n data matrix
#' @param matJointV n x sum(r_t) joint structure basis matrix
#' @param matJointOrder List that keeps the order of joint blocks in matJointV
#' @param matJointRanks Vector containing the rank of each joint block in matJointV
#'
#' @return A list containing `matBlockMap` (joint block matrix map) and `matLoadingMap` (joint block loading map)
MatReconstructMJ <- function(X, matJointV, matJointOrder, matJointRanks) {
  # A data block can legitimately have no retained joint structure.  Return
  # empty component maps instead of passing a NULL/zero-column basis to QR.
  if (is.null(matJointV) || ncol(matJointV) == 0L || length(matJointRanks) == 0L) {
    return(list(matBlockMap = list(), matLoadingMap = list()))
  }

  # Solve linear equation to get the loadings
  #matLoading <- t(solve(t(matJointV), t(X)))
  #matLoading <- t(solve(t(matJointV) %*% matJointV) %*% t(matJointV) %*% t(X))
  #matLoading <- t(ginv(matJointV) %*% t(X))
  matLoading <- t(qr.solve(matJointV, t(X)))

  #print(matLoading)
  matBlockMap <- list()
  matLoadingMap <- list()
  cum <- 0


  for (j in seq_along(matJointOrder)) {
    t <- matJointOrder[[j]]
    r <- matJointRanks[j]
    Vhat <- matJointV[, (cum + 1):(cum + r)]
    #Bhat <- matLoading[, (cum + 1):(cum + r)]
    Bhat <- matLoading[, (cum + 1):(cum + r)]
    matBlockMap[[as.character(t)]] <- Bhat %*% t(Vhat)
    matLoadingMap[[as.character(t)]] <- Bhat
    cum <- cum + r
  }

  return(list(matBlockMap = matBlockMap, matLoadingMap = matLoadingMap))
}







#' DJIVEReconstructMJ - Reconstruct joint blocks from data blocks
#'
#' This function reconstructs joint blocks from the provided data matrices
#' using the estimated partially shared structures.
#'
#' @param datablock List of d_k x n data matrices.
#' @param dataname Vector of data matrix names.
#' @param outMap List mapping between block index set and estimated partially shared structure.
#' @param keyIdxMap List mapping between joint block index and data blocks.
#' @param jointBlockOrder Vector to record the joint block order.
#' @param doubleCenter Flag indicating whether to perform double centering.
#' @param categorical_blocks Optional list aligned with datablock; entries are
#'   observed categorical blocks (raw or tagged), with NULL for continuous blocks.
#'
#' @return A list containing the joint basis map, loadings, and joint blocks.
#' @export
DJIVEReconstructMJ <- function(datablock, dataname, outMap, keyIdxMap, jointBlockOrder, doubleCenter, categorical_blocks = NULL) {
  nb <- length(datablock)
  if (is.null(categorical_blocks)) categorical_blocks <- rep(list(NULL), nb)
  if (length(categorical_blocks) != nb) stop("categorical_blocks must align with datablock.")
  categoricalReconstruction <- setNames(vector("list", nb), dataname)
  matJointV <- vector("list", nb)
  matJointOrder <- vector("list", nb)
  matJointRanks <- vector("list", nb)

  # # Initialize joint basis, order, and ranks
  # for (ib in seq_len(nb)) {
  #   matJointV[[ib]] <- matrix()
  #   matJointOrder[[ib]] <- character()
  #   matJointRanks[[ib]] <- numeric()
  # }

  for (ib in seq_len(nb)) {
    # Joint bases live in sample space, so they have one row per sample.
    matJointV[[ib]] <- matrix(0, nrow = ncol(datablock[[ib]]), ncol = 0)
    matJointOrder[[ib]] <- character()
    matJointRanks[[ib]] <- numeric()
  }

  # print(matJointV)
  #print(dim(matJointV))

  # Re-rotate joint bases for interpretation
  commonNormalizedBasisMap <- list()
  for (idx in 1:(2^nb - 1)) {
    idxStr <- as.character(idx)
    if (idxStr %in% names(outMap)) {
      blockIn <- keyIdxMap[[idxStr]]
      jointSpace <- outMap[[idxStr]]
      jointSpaceProj <- jointSpace %*% t(jointSpace)
      jointDim <- ncol(jointSpace)
      dataStack <- do.call(rbind, datablock[blockIn])
      projected <- dataStack %*% jointSpaceProj
      commonNormalizedBasisMap[[idxStr]] <- if (min(dim(projected)) < 3L ||
                                               jointDim >= min(dim(projected))) {
        svd(projected, nu = 0L, nv = jointDim)$v
      } else svds(projected, jointDim)$v
    }
  }
  outMapRaw <- outMap
  outMap <- commonNormalizedBasisMap
  #print(dim(outMap))
  # Collect joint block basis for each data matrix
  for (i in seq_along(jointBlockOrder)) {
    t <- jointBlockOrder[[i]]
    V <- outMap[[t]]

    # print(V)
    # print(dim(V))
    #print(matJointV)
    #print(dim(matJointV))

    r <- ncol(V)
    blockIdx <- keyIdxMap[[t]]
    cat(sprintf("The rank of joint block among %s: %d\n", paste(dataname[blockIdx], collapse = ", "), r))
    for (ib in blockIdx) {
      matJointV[[ib]] <- cbind(matJointV[[ib]], V)
      matJointOrder[[ib]] <- c(matJointOrder[[ib]], t)
      matJointRanks[[ib]] <- c(matJointRanks[[ib]], r)
    }
  }

  # Calculate the loadings and the joint block in each data matrix
  matLoadings <- vector("list", nb)
  matBlocks <- vector("list", nb)
  for (ib in seq_len(nb)) {
    datablockc <- datablock[[ib]]
    if (doubleCenter == 1 && is.null(categorical_blocks[[ib]])) {
      bsize <- dim(datablock[[ib]])
      d <- bsize[1]
      n <- bsize[2]
      datablockc <- datablock[[ib]] - rowMeans(datablock[[ib]]) - matrix(colMeans(datablock[[ib]]), nrow = d, ncol = n, byrow = TRUE) + mean(datablock[[ib]])
    }
    res <- if (is.null(categorical_blocks[[ib]])) {
      MatReconstructMJ(datablockc, matJointV[[ib]], matJointOrder[[ib]], matJointRanks[[ib]])
    } else {
      choir_reconstruct(categorical_blocks[[ib]], matJointV[[ib]],
                        matJointOrder[[ib]], matJointRanks[[ib]])
    }
    if (!is.null(categorical_blocks[[ib]])) categoricalReconstruction[[ib]] <- res
    matBlocks[[ib]] <- res$matBlockMap
    matLoadings[[ib]] <- res$matLoadingMap
  }

  # Output structure
  outstruct <- list(
    jointBasisMap = outMap,
    matLoadings = matLoadings,
    matBlocks = matBlocks,
    keyIdxMap = keyIdxMap,
    categoricalReconstruction = categoricalReconstruction
  )
  return(outstruct)
}
