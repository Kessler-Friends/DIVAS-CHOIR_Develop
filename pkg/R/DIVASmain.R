#' Data integration via analysis of subspaces
#'
#' Main function for DIVAS analysis. Given a list of data blocks with matched columns (samples), will return identified joint
#' structure with diagnostic plots.
#'
#' @param continuous_datablock A list of matrices with the same number of columns (samples).
#' @param categorical_datablock A list of categorical matrices or data.frames. For data.frames, each column is treated as a categorical variable measured across samples and is one-hot encoded before analysis.
#' @param nsim Number of bootstrap resamples for inferring angle bounds.
#' @param iprint Whether to print diagnostic figures.
#' @param colCent Whether to column centre the input data blocks.
#' @param rowCent Whether to row centre the input data blocks.
#' @param figdir If not NULL, then diagnostic plots will be saved to this directory.
#' @param seed Optional. An integer to set the seed for the random number generator to ensure reproducibility of the bootstrap analysis. Default is `NULL`.
#'
#' @return A list containing DIVAS integration results. Most important ones include
#'   \describe{
#'     \item{matBlocks}{List of scores representing shared and partially shared joint structures.}
#'     \item{matLoadings}{List of loadings linking features in each data block with scores.}
#'     \item{keyIdxMap}{Mapping between indices of the previous lists and data blocks.}
#'   }
#'   See Details for more explanations.
#'
#' @details
#' DIVASmain returns a list containing all important information returned from the DIVAS algorithm.
#' For users, the most important ones are scores (matBlocks), loadings (matLoadings) and an index
#' book (keyIdxMap) explaining what joint structures each score or loading matrix correpsond to.
#'
#' matBlocks is a list containing scores. Each element of matBlocks is indexed by a number.
#' For example, suppose one of the indices is "7", then keyIdxMap[["7"]] contains indices of data blocks
#' corresponding to the index 7. That is, matBlocks[["7"]] contains the scores for all samples
#' representing the joint structures of data blocks in keyIdxMap[["7"]].
#'
#' @references
#' Prothero, J., Jiang, M., Hannig, J., Tran-Dinh, Q., Ackerman, A. and Marron, J. S. (2024).
#' Data integration via analysis of subspaces (DIVAS). Test.
#'
DIVASmain <- function(
    continuous_datablock = NULL, categorical_datablock = NULL, nsim = 400,
    iprint = TRUE, colCent = FALSE, rowCent = FALSE, figdir = NULL, seed = NULL
  ){
  continuous_datablock <- if (is.null(continuous_datablock)) list() else continuous_datablock
  categorical_datablock <- if (is.null(categorical_datablock)) list() else categorical_datablock

  # Initialize continuous parameters
  continuous_nb <- length(continuous_datablock)
  continuous_dataname <- names(continuous_datablock)
  if(is.null(continuous_dataname) && continuous_nb > 0){
    warning("Input continuous_datablock is unnamed, generic names for continuous data blocks generated.")
    continuous_dataname <- paste0("ContinuousDatablock", 1:continuous_nb)
  }

  # Initialize categorical parameters
  categorical_nb <- length(categorical_datablock)
  categorical_dataname <- names(categorical_datablock)
  if(is.null(categorical_dataname) && categorical_nb > 0){
    warning("Input categorical_datablock is unnamed, generic names for categorical data blocks generated.")
    categorical_dataname <- paste0("CategoricalDatablock", 1:categorical_nb)
  }

  # Some tuning parameters for algorithms
  theta0 <- 45
  optArgin <- list(0.5, 1000, 1.05, 50, 1e-3, 1e-3)
  filterPerc <- 1 - (2 / (1 + sqrt(5))) # "Golden Ratio"
  continuous_noisepercentile <- rep(0.5, continuous_nb)


  continuous_rowSpaces <- vector("list", continuous_nb)
  # datablockc <- vector("list", nb)
  for (ib in seq_len(continuous_nb)) {
    continuous_rowSpaces[[ib]] <- 0
    continuous_datablock[[ib]] <- MatCenterJP(continuous_datablock[[ib]], colCent, rowCent)
  }

  categorical_rowSpaces <- vector("list", categorical_nb)
  for (ib in seq_len(categorical_nb)) {
    categorical_rowSpaces[[ib]] <- 0
    categorical_datablock[[ib]] <- MatCenterDIVASCHOIR(categorical_datablock[[ib]])
  }

  # Step 1: Estimate continuous signal space and perturbation angle
  Phase1_continuous <- DJIVESignalExtractJP(
    datablock = continuous_datablock, nsim = nsim,
    iplot = FALSE, colCent = colCent, rowCent = rowCent, cull = filterPerc, noisepercentile = continuous_noisepercentile,
    seed = seed
  )

  # Step 1: Estimate categorical signal space
  if (categorical_nb > 0) {
    Phase1_categorical <- DIVASCHOIRSignalExtract(
      datablock = categorical_datablock,
      nsim = nsim,
      iplot = FALSE,
      cull = filterPerc,
      seed = seed
    )
  } else {
    Phase1_categorical <- list(
      VBars = list(), UBars = list(), phiBars = numeric(), psiBars = numeric(),
      EHats = list(), rBars = numeric(), singVals = list(), singValsHat = list(),
      rSteps = list(), VVHatCacheBars = list(), UUHatCacheBars = list()
    )
  }

  Phase1_combined <- list(
    VBars = c(Phase1_continuous$VBars, Phase1_categorical$VBars),
    UBars = c(Phase1_continuous$UBars, Phase1_categorical$UBars),
    phiBars = c(Phase1_continuous$phiBars, Phase1_categorical$phiBars),
    psiBars = c(Phase1_continuous$psiBars, Phase1_categorical$psiBars),
    EHats = c(Phase1_continuous$EHats, Phase1_categorical$EHats),
    rBars = c(Phase1_continuous$rBars, Phase1_categorical$rBars),
    singVals = c(Phase1_continuous$singVals, Phase1_categorical$singVals),
    singValsHat = c(Phase1_continuous$singValsHat, Phase1_categorical$singValsHat),
    rSteps = c(Phase1_continuous$rSteps, Phase1_categorical$rSteps),
    VVHatCacheBars = c(Phase1_continuous$VVHatCacheBars, Phase1_categorical$VVHatCacheBars),
    UUHatCacheBars = c(Phase1_continuous$UUHatCacheBars, Phase1_categorical$UUHatCacheBars)
  )

  combined_datablock <- c(continuous_datablock, categorical_datablock)
  combined_dataname <- c(continuous_dataname, categorical_dataname)

  # VBars <- Phase1[[1]]
  # UBars <- Phase1[[2]]
  # phiBars <- Phase1[[3]]
  # psiBars <- Phase1[[4]]
  # rBars <- Phase1[[6]]
  # VVHatCacheBars <- Phase1[[10]]
  # UUHatCacheBars <- Phase1[[11]]


  # Step 2: Estimate joint and partially joint structure
  Phase2 <- DJIVEJointStrucEstimateJP(
    VBars = Phase1_combined$VBars, UBars = Phase1_combined$UBars, phiBars =  Phase1_combined$phiBars, psiBars =  Phase1_combined$psiBars,
    rBars = Phase1_combined$rBars, dataname = combined_dataname, theta0 = theta0,
    optArgin = optArgin, iprint = iprint, figdir = figdir
  )

  # outMap <- Phase2[[1]]
  # keyIdxMap <- Phase2[[2]]
  # jointBlockOrder <- Phase2[[4]]

  # Step 3: Reconstruct DJIVE decomposition
  outstruct <- DJIVEReconstructMJ(
    datablock = combined_datablock, dataname =  combined_dataname, outMap =  Phase2$outMap,
    keyIdxMap =  Phase2$keyIdxMap, jointBlockOrder =  Phase2$jointBlockOrder, doubleCenter =  0
  )

  outstruct$rBars <- Phase1_combined$rBars
  outstruct$phiBars <- Phase1_combined$phiBars
  outstruct$psiBars <- Phase1_combined$psiBars
  outstruct$VBars <- Phase1_combined$VBars
  outstruct$UBars <- Phase1_combined$UBars
  outstruct$VVHatCacheBars <- Phase1_combined$VVHatCacheBars
  outstruct$UUHatCacheBars <- Phase1_combined$UUHatCacheBars
  outstruct$Phase1_continuous <- Phase1_continuous
  outstruct$Phase1_categorical <- Phase1_categorical
  outstruct$Phase1_combined <- Phase1_combined
  outstruct$jointBasisMapRaw <- Phase2$outMap

  # Automatically generate keymapname from keymapid
  ids <- as.integer(names(outstruct$keyIdxMap))
  num_blocks <- length(combined_dataname)
  
  keymapname <- sapply(ids, function(id) {
    binary_str <- R.utils::intToBin(id)
    padded_binary_str <- sprintf(paste0("%0", num_blocks, "s"), binary_str)
    binary_chars <- strsplit(padded_binary_str, "")[[1]]
    selected_indices <- which(rev(binary_chars) == '1')
    selected_names <- combined_dataname[selected_indices]
    paste(selected_names, collapse = "+")
  })
  
  names(keymapname) <- names(outstruct$keyIdxMap)
  outstruct$keymapname <- keymapname
  outstruct$combined_datablock <- combined_datablock
  outstruct$combined_dataname <- combined_dataname
  outstruct$categorical_datablock <- categorical_datablock
  outstruct$categorical_dataname <- categorical_dataname

  diagnostic_plots <- tryCatch(
    DJIVEAngleDiagnosticJP(
      datablock = combined_datablock,
      dataname = combined_dataname,
      outstruct = outstruct,
      randseed = if (is.null(seed)) 1 else seed,
      titlestr = "DIVAS diagnostics"
    ),
    error = function(e) {
      warning("Failed to generate DIVAS diagnostic plots: ", conditionMessage(e))
      NULL
    }
  )

  outstruct$diagnostic_plots <- diagnostic_plots
  if (!is.null(diagnostic_plots) && !is.null(figdir)) {
    dir.create(figdir, recursive = TRUE, showWarnings = FALSE)
    ggplot2::ggsave(
      filename = file.path(figdir, "rank_breakdown.png"),
      plot = diagnostic_plots$rank,
      width = 12,
      height = 6,
      dpi = 300
    )
    ggplot2::ggsave(
      filename = file.path(figdir, "joint_structure_score_diagnostics.png"),
      plot = diagnostic_plots$score,
      width = 12,
      height = 8,
      dpi = 300
    )
    ggplot2::ggsave(
      filename = file.path(figdir, "joint_structure_loading_diagnostics.png"),
      plot = diagnostic_plots$loading,
      width = 12,
      height = 8,
      dpi = 300
    )
  }

  return(outstruct)
}
