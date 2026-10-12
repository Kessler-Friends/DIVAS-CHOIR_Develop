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
#' @param ReturnDetail Logical. If FALSE (default), return Scores, scoresList, Loadings, and LoadingsNames. If TRUE, also return the detailed decomposition and diagnostic plots.
#' @param datablock Optional legacy alias for continuous_datablock.
#' @param categorical_rank Optional initial rank for categorical blocks. The
#'   original DIVAS angle rule then culls unstable directions by default.
#' @param categorical_method `"fixed_rank"` (default) or `"nuclear_cv"`. The latter selects a nuclear-norm penalty by cross-validation of held-out categorical responses.
#' @param categorical_cv_folds Number of folds for categorical nuclear-norm cross-validation.
#' @param categorical_lambda_grid Optional positive penalty values to consider. By default, ten values are generated from the MCA singular spectrum.
#' @param categorical_cull Apply the original DIVAS random-angle rank culling rule to categorical directions. Defaults to TRUE.
#'
#' @return A list containing DIVAS integration results. Most important ones include
#'   \describe{
#'     \item{matBlocks}{List of scores representing shared and partially shared joint structures.}
#'     \item{matLoadings}{List of loadings linking features in each data block with scores.}
#'     \item{keyIdxMap}{Mapping between indices of the previous lists and data blocks.}
#'   }
#'   The compact result contains Scores (combined sample bases), scoresList
#'   (recalculated scores per block), Loadings (normalized inverse loadings),
#'   and LoadingsNames (component names per block). Set ReturnDetail = TRUE
#'   to include the internal results and diagnostic_plots.
#'   See Details for more explanations.
#'
#' @details
#' With ReturnDetail = TRUE, DIVASmain also returns the internal results from the DIVAS algorithm.
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
    iprint = TRUE, colCent = FALSE, rowCent = FALSE, figdir = NULL, seed = NULL,
    categorical_rank = NULL, categorical_method = c("fixed_rank", "nuclear_cv"),
    categorical_cv_folds = 5L, categorical_lambda_grid = NULL,
    categorical_cull = TRUE, ReturnDetail = FALSE, datablock = NULL
  ){
  if (!is.null(datablock)) {
    if (!is.null(continuous_datablock)) {
      stop("Supply only one of `datablock` and `continuous_datablock`.")
    }
    continuous_datablock <- datablock
  }
  categorical_method <- match.arg(categorical_method)
  if (categorical_method == "nuclear_cv" && !is.null(categorical_rank)) {
    stop("`categorical_rank` cannot be combined with `categorical_method = 'nuclear_cv'`.")
  }
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
  categorical_cv <- vector("list", categorical_nb)
  categorical_lambda <- rep(NA_real_, categorical_nb)
  for (ib in seq_len(categorical_nb)) {
    categorical_rowSpaces[[ib]] <- 0
    if (categorical_method == "nuclear_cv") {
      categorical_cv[[ib]] <- categorical_nuclear_cv(
        categorical_datablock[[ib]], folds = categorical_cv_folds,
        lambda_grid = categorical_lambda_grid,
        seed = if (is.null(seed)) NULL else seed + ib - 1L
      )
      categorical_lambda[ib] <- categorical_cv[[ib]]$lambda
    }
    categorical_datablock[[ib]] <- MatCenterDIVASCHOIR(categorical_datablock[[ib]])
  }

  # Step 1: Estimate continuous signal space and perturbation angle
  if (continuous_nb > 0) {
    Phase1_continuous <- DJIVESignalExtractJP(
      datablock = continuous_datablock, nsim = nsim,
      iplot = FALSE, colCent = colCent, rowCent = rowCent, cull = filterPerc, noisepercentile = continuous_noisepercentile,
      seed = seed
    )
  } else {
    Phase1_continuous <- list(
      VBars = list(), UBars = list(), phiBars = numeric(), psiBars = numeric(),
      EHats = list(), rBars = numeric(), singVals = list(), singValsHat = list(),
      rSteps = list(), VVHatCacheBars = list(), UUHatCacheBars = list()
    )
  }

  # Step 1: Estimate categorical signal space
  if (categorical_nb > 0) {
    Phase1_categorical <- DIVASCHOIRSignalExtract(
      datablock = categorical_datablock,
      nsim = nsim,
      iplot = FALSE,
      cull = filterPerc,
      seed = seed,
      fixed_rank = categorical_rank,
      lambda = if (categorical_method == "nuclear_cv") categorical_lambda else NULL,
      cull_rank = categorical_cull
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

  # Original post-processing operates on the combined feature-by-sample blocks.
  datablock <- combined_datablock
  dataname <- combined_dataname
  nb <- length(datablock)
  sample_ids <- if (nb > 0) colnames(datablock[[1]]) else NULL

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
  outstruct$categorical_method <- categorical_method
  outstruct$categorical_lambda <- if (categorical_method == "nuclear_cv") categorical_lambda else NULL
  outstruct$categorical_cv <- if (categorical_method == "nuclear_cv") categorical_cv else NULL
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

  # ---- Attach sample IDs to results ----
  if (!is.null(sample_ids)) {
    if (iprint) {
      cat("Attaching sample IDs to score matrices (jointBasisMap and indivBasisMap)...\n")
    }
    # For joint structures
    if (length(outstruct$jointBasisMap) > 0) {
      outstruct$jointBasisMap <- lapply(outstruct$jointBasisMap, function(mat) {
        if (!is.null(mat)) rownames(mat) <- sample_ids
        return(mat)
      })
    }
    # For individual structures
    if (length(outstruct$indivBasisMap) > 0) {
      outstruct$indivBasisMap <- lapply(outstruct$indivBasisMap, function(mat) {
        if (!is.null(mat)) rownames(mat) <- sample_ids
        return(mat)
      })
    }
  }

  # ---- Create comprehensive sample-score matrix ----
  if (iprint) {
    cat("Creating comprehensive sample-score matrix...\n")
  }

  # Initialize list to store all score columns
  score_columns <- list()
  column_names <- c()

  # Process joint structures (sorted by number of blocks, descending)
  if (length(outstruct$jointBasisMap) > 0) {
    # Get joint structure information
    joint_info <- data.frame(
      id = names(outstruct$jointBasisMap),
      stringsAsFactors = FALSE
    )

    # Calculate number of blocks for each joint structure
    joint_info$num_blocks <- sapply(joint_info$id, function(id) {
      length(outstruct$keyIdxMap[[id]])
    })

    # Sort by number of blocks (descending) then by id
    joint_info <- joint_info[order(-joint_info$num_blocks, joint_info$id), ]

    # Process each joint structure
    for (i in seq_len(nrow(joint_info))) {
      id <- joint_info$id[i]
      num_blocks <- joint_info$num_blocks[i]
      mat <- outstruct$jointBasisMap[[id]]

      if (!is.null(mat) && ncol(mat) > 0) {
        # Determine if this is a true joint structure (>1 block) or individual structure
        if (num_blocks > 1) {
          # True joint structure: use actual block names like "RNA+PRO+MIC-1", "RNA+PRO-1", etc.
          block_indices <- outstruct$keyIdxMap[[id]]
          block_names <- dataname[block_indices]
          joint_name <- paste(block_names, collapse = "+")

          for (rank_idx in seq_len(ncol(mat))) {
            col_name <- paste0(joint_name, "-", rank_idx)
            column_names <- c(column_names, col_name)
            score_columns[[col_name]] <- mat[, rank_idx]
          }
        } else {
          # Individual structure stored in jointBasisMap: "RNA-Individual-1", etc.
          block_idx <- outstruct$keyIdxMap[[id]][1]
          block_name <- dataname[block_idx]

          for (rank_idx in seq_len(ncol(mat))) {
            col_name <- paste0(block_name, "-Individual-", rank_idx)
            column_names <- c(column_names, col_name)
            score_columns[[col_name]] <- mat[, rank_idx]
          }
        }
      }
    }
  }

  # Process individual structures (if they exist separately)
  if (length(outstruct$indivBasisMap) > 0) {
    for (id in names(outstruct$indivBasisMap)) {
      mat <- outstruct$indivBasisMap[[id]]

      if (!is.null(mat) && ncol(mat) > 0) {
        # Get block name for individual structure
        block_idx <- outstruct$keyIdxMap[[id]][1]  # Individual structures have only one block
        block_name <- dataname[block_idx]

        # Create column names like "RNA-Individual-1", "PRO-Individual-2", etc.
        for (rank_idx in seq_len(ncol(mat))) {
          col_name <- paste0(block_name, "-Individual-", rank_idx)
          column_names <- c(column_names, col_name)
          score_columns[[col_name]] <- mat[, rank_idx]
        }
      }
    }
  }

  # Create the comprehensive matrix
  if (length(score_columns) > 0) {
    # Combine all score columns into a matrix
    sample_score_matrix <- do.call(cbind, score_columns)

    # Debug information
    if (iprint) {
      cat("Debug: score_columns length:", length(score_columns), "\n")
      cat("Debug: column_names length:", length(column_names), "\n")
      cat("Debug: matrix dimensions:", dim(sample_score_matrix), "\n")
    }

    # Ensure column names match matrix dimensions
    if (ncol(sample_score_matrix) == length(column_names)) {
      colnames(sample_score_matrix) <- column_names
    } else {
      if (iprint) {
        cat("Warning: Column names length (", length(column_names),
            ") does not match matrix columns (", ncol(sample_score_matrix), ")\n")
      }
      # Generate generic column names if there's a mismatch
      colnames(sample_score_matrix) <- paste0("Component_", seq_len(ncol(sample_score_matrix)))
    }

    # Set row names to sample IDs if available
    if (!is.null(sample_ids)) {
      rownames(sample_score_matrix) <- sample_ids
    }

    # Add to output structure
    outstruct$sampleScoreMatrix <- sample_score_matrix

    if (iprint) {
      cat("Sample-score matrix created with dimensions:", dim(sample_score_matrix), "\n")
      cat("Columns:", paste(colnames(sample_score_matrix), collapse = ", "), "\n")
    }
  } else {
    outstruct$sampleScoreMatrix <- NULL
    if (iprint) {
      cat("No score matrices found to create sample-score matrix.\n")
    }
  }

  # ---- Compute inverse loadings and recalculated scores ----
  if (iprint) {
    cat("Computing inverse loadings and recalculated scores...\n")
  }

  # Initialize output structures
  outstruct$Loadings <- vector("list", nb)
  outstruct$scoresList <- vector("list", nb)
  outstruct$LoadingsNames <- vector("list", nb)
  outstruct$names_vec <- vector("list", nb)
  names(outstruct$Loadings) <- dataname
  names(outstruct$scoresList) <- dataname
  names(outstruct$LoadingsNames) <- dataname
  names(outstruct$names_vec) <- dataname

  # Helper function: computLoadings (exactly from CaseCOVID.Rmd)
  computLoadings <- function(divasRes, omic_block = 1) {

    # Step 1: Concatenate all loadings for the specified omic block
    all_loadings_list <- divasRes$matLoadings[[omic_block]]
    # Remove any NULL components
    valid_loadings <- Filter(function(mat) !is.null(mat) && ncol(as.matrix(mat)) > 0,
                             all_loadings_list)
    newNames <- names(valid_loadings)
    newNames <- divasRes$keymapname[newNames]
    names(valid_loadings) <- newNames

    if (length(valid_loadings) == 0) {
      stop("No valid loadings found for this omic block")
    }


    L_concatenated <- do.call(cbind, valid_loadings)

    # Step 3: Compute generalized inverse of concatenated loadings
    L_concat_pinv_T <- t(MASS::ginv(L_concatenated))
    loading_norms <- sqrt(colSums(L_concat_pinv_T^2))
    loading_norms[loading_norms == 0] <- 1
    L_concat_pinv_T <- scale(L_concat_pinv_T, center = FALSE, scale = loading_norms)


    # Step 4: Add rownames and column names
    rownames(L_concat_pinv_T) <- rownames(L_concatenated)

    # Create column names based on data block combinations and component indices
    col_names <- c()
    for (i in seq_along(valid_loadings)) {
      loading_name <- names(valid_loadings)[i]  # This is the keymapname (e.g., "RNA+Protein")
      n_components <- ncol(as.matrix(valid_loadings[[i]]))

      # Generate column names like "RNA+Protein-1", "RNA+Protein-2", etc.
      component_names <- paste0(loading_name, "-", seq_len(n_components))
      col_names <- c(col_names, component_names)
    }
    # Set column names
    colnames(L_concat_pinv_T) <- col_names

    return(L_concat_pinv_T)
  }

  # Process each data block using computLoadings
  for (block_idx in seq_len(nb)) {
    if (iprint) {
      cat("Processing block", block_idx, ":", dataname[block_idx], "...\n")
    }

    # MASS is imported, so we can use it directly

    # Use computLoadings function
    tryCatch({
      newLoadings <- computLoadings(outstruct, omic_block = block_idx)

      # Store the inverse loadings as a matrix per block:

      # - More intuitive structure than nested lists
      outstruct$Loadings[[block_idx]] <- newLoadings
      outstruct$LoadingsNames[[block_idx]] <- colnames(newLoadings)
      outstruct$names_vec[[block_idx]] <- colnames(newLoadings)

      # Compute recalculated scores following CaseCOVID.Rmd:
      # Shat <- crossprod(XtrainCent, newLoadings)
      X_block <- as.matrix(datablock[[block_idx]])
      if (!is.null(sample_ids)) {
        X_block <- X_block[, sample_ids, drop = FALSE]
      }

      # Categorical blocks already use the centered MCA representation.
      X_block_cent <- if (block_idx <= continuous_nb) {
        X_block - rowMeans(X_block)
      } else {
        X_block
      }
      S_recalc <- crossprod(X_block_cent, newLoadings)

      score_norms <- sqrt(colSums(S_recalc^2))
      score_norms[score_norms == 0] <- 1
      S_recalc <- sweep(S_recalc, 2, score_norms, "/")

      # Add sample names as rownames and column names
      if (!is.null(sample_ids)) {
        rownames(S_recalc) <- sample_ids
      }
      colnames(S_recalc) <- colnames(newLoadings)

      # Store the recalculated scores
      outstruct$scoresList[[block_idx]] <- S_recalc

      if (iprint) {
        cat("  Inverse loadings computed with dimensions:", dim(newLoadings), "\n")
        cat("  Recalculated scores computed with dimensions:", dim(S_recalc), "\n")
      }
    }, error = function(e) {
      if (iprint) {
        cat("  Error processing block", block_idx, ":", e$message, "\n")
      }
      outstruct$Loadings[block_idx] <<- list(NULL)
      outstruct$scoresList[block_idx] <<- list(NULL)
    })
  }

  if (iprint) {
    cat("DIVAS is complete.\n")
  }


  if (!ReturnDetail) {
    # Compact result with updated naming and order per user requirements
    compact <- list(
      Scores = outstruct$sampleScoreMatrix,  # first: rename sampleScoreMatrix to Scores
      scoresList = outstruct$scoresList,     # second: directly use scoresList (no renaming needed)
      Loadings = outstruct$Loadings,         # third: keep Loadings as matrix format
      LoadingsNames = outstruct$LoadingsNames # fourth: keep LoadingsNames (remove names_vec)
    )
    return(compact)
  } else {
    # Full details with updated naming and priority order
    prioritized <- list(
      Scores = outstruct$sampleScoreMatrix,  # first: rename sampleScoreMatrix to Scores
      scoresList = outstruct$scoresList,     # second: directly use scoresList (no renaming needed)
      Loadings = outstruct$Loadings,         # third: keep Loadings as matrix format
      LoadingsNames = outstruct$LoadingsNames # fourth: keep LoadingsNames
    )
    # Append the rest of the fields preserving originals (excluding renamed and duplicate ones)
    remaining_names <- setdiff(names(outstruct), c("sampleScoreMatrix", "scoresList", "Loadings", "LoadingsNames", "names_vec", "jointBasisMap"))
    for (nm in remaining_names) {
      prioritized[[nm]] <- outstruct[[nm]]
    }
    return(prioritized)
  }
}
