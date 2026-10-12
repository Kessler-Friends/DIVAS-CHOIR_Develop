#!/usr/bin/env Rscript
# Convert X3 of the original DIVAS toy data into 36 binary variables.
# X1 and X2 remain the original continuous matrices. The conversion uses X3
# alone, so the other blocks cannot influence the categorical measurements.
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Run with Rscript")
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg)))
pkg_dir <- dirname(script_dir)
lib <- file.path(Sys.getenv("LOCALAPPDATA"), "R", "mfa_methods_library")
if (dir.exists(lib)) .libPaths(c(lib, .libPaths()))
if (!requireNamespace("R.matlab", quietly = TRUE)) stop("Install R.matlab")

source_file <- file.path(pkg_dir, "inst", "extdata", "toyDataThreeWay.mat")
toy <- suppressWarnings(R.matlab::readMat(source_file))
blocks <- lapply(seq_len(3L), function(i) as.matrix(toy$datablock[[i]][[1L]]))
if (!identical(lapply(blocks, dim), list(c(200L, 400L), c(400L, 400L),
                                         c(10000L, 400L))))
  stop("Unexpected original toy block dimensions")

n <- ncol(blocks[[3L]])
ids <- sprintf("S%03d", seq_len(n))
for (i in seq_along(blocks)) colnames(blocks[[i]]) <- ids
rownames(blocks[[1L]]) <- sprintf("X1_%03d", seq_len(nrow(blocks[[1L]])))
rownames(blocks[[2L]]) <- sprintf("X2_%03d", seq_len(nrow(blocks[[2L]])))

# The original X3 has three separated leading singular values. Its leading
# sample scores define three axes. Twelve slightly perturbed projections of
# each axis give balanced categorical measurements of all three directions.
# No other block is used in this transformation.
eig <- eigen(crossprod(blocks[[3L]]), symmetric = TRUE)
source_scores <- eig$vectors[, seq_len(3L), drop = FALSE]
for (j in seq_len(3L)) {
  anchor <- which.max(abs(source_scores[, j]))
  if (source_scores[anchor, j] < 0) source_scores[, j] <- -source_scores[, j]
}
source_scores <- source_scores * sqrt(n)
rownames(source_scores) <- ids
colnames(source_scores) <- paste0("X3_PC", seq_len(3L))

set.seed(1207L)
projection_vectors <- do.call(rbind, lapply(seq_len(3L), function(axis) {
  base <- matrix(rep(diag(3L)[axis, ], each = 12L), nrow = 12L)
  base + matrix(stats::rnorm(36L, sd = 0.2), nrow = 12L)
}))
projection_vectors <- projection_vectors /
  sqrt(rowSums(projection_vectors^2))
rownames(projection_vectors) <- sprintf("Trait_%02d", seq_len(nrow(projection_vectors)))
colnames(projection_vectors) <- colnames(source_scores)
projections <- source_scores %*% t(projection_vectors)
categorical <- as.data.frame(lapply(seq_len(ncol(projections)), function(j) {
  factor(ifelse(projections[, j] <= stats::median(projections[, j]),
                "Low", "High"), levels = c("Low", "High"))
}), stringsAsFactors = TRUE)
names(categorical) <- rownames(projection_vectors)
rownames(categorical) <- ids
if (anyNA(categorical) || any(vapply(categorical, nlevels, integer(1)) != 2L))
  stop("Categorical conversion failed")

dataset <- list(
  continuous = list(X1 = blocks[[1L]], X2 = blocks[[2L]]),
  categorical = list(X3_categories = categorical),
  source_x3_scores = source_scores,
  source_x3_singular_values = sqrt(pmax(eig$values, 0)),
  projection_vectors = projection_vectors,
  source_file = "toyDataThreeWay.mat",
  conversion = "Three X3 sample PCs, 12 fixed projections per PC, median split into Low/High"
)
output_file <- file.path(pkg_dir, "inst", "extdata", "toyDataThreeWayCategorical.rds")
saveRDS(dataset, output_file, compress = "xz")
message("Saved mixed toy data: ", output_file)
