\name{DIVASCHOIRSignalExtract}
\alias{DIVASCHOIRSignalExtract}
\title{Categorical signal extraction for DIVAS-CHOIR}
\usage{
DIVASCHOIRSignalExtract(datablock, nsim = 400, iplot = FALSE,
  cull = 0.5, seed = NULL, fixed_rank = NULL,
  cull_rank = FALSE, method = c("mca", "nuclear_norm"), lambda_grid = NULL,
  simulate_null = TRUE, bootstrap_tol = 1e-7, null_nsim = 1000L)
}
\arguments{
\item{simulate_null}{Compute null-angle diagnostics; required when rank culling is enabled.}
\item{bootstrap_tol}{Positive nuclear refit convergence tolerance. Original/CV fits use 1e-7.}
\item{null_nsim}{Positive integer number of categorical null simulations, independent of nsim bootstrap samples. Defaults to 1000.}
\item{datablock}{A list of raw categorical data frames (subjects by variables), or tagged matrices prepared for the selected method. The extractor prepares original and bootstrap inputs.}
\item{nsim}{Number of bootstrap samples.}
\item{iplot}{Whether to plot singular values.}
\item{cull}{Multiplier applied to the fifth percentile of the simulated categorical null angles for rank culling. Defaults to 0.5.}
\item{seed}{Optional random seed.}
\item{fixed_rank}{Optional initial latent signal rank, either shared or one value per block, before culling. Nuclear fits are truncated to the leading singular components at this rank.}
\item{method}{Estimator used for observed, null and bootstrap data: "mca" or "nuclear_norm". Nuclear norm fits multinomial likelihood with a nuclear-norm penalty on centered eta. Five-fold response CV selects lambda and rank. The selected lambda stays fixed in null and bootstrap fits; their fitted ranks may differ.}
\item{lambda_grid}{Optional finite nonnegative penalties for nuclear five-fold response CV selection. NULL uses factors 1.05, 0.75, 0.45, 0.25, and 0.1 times the intercept-only score spectral norm. Ignored for MCA.}
\item{cull_rank}{Apply rank culling using simulated categorical null angles. Defaults to FALSE, retaining all initial signal directions. No direction is forced to remain. If none pass, return rank zero, empty bases and caches, and 90-degree angle placeholders.}
}
\value{
A list containing categorical signal bases, ranks, singular values, bootstrap angle bounds, and \code{nullAngleCaches}, plus \code{estimation} with method, rank and penalty diagnostics. Each null cache contains null_nsim sample-space (\code{subject}) and category-space (\code{feature}) angles in degrees when calibration is performed.
}
\description{
Computes categorical signal spaces and bootstrap angle bounds with method-specific preparation: centered/scaled indicators for MCA, unchanged one-hot data and latent centering for nuclear likelihood.
}
\details{
Nuclear selection uses five-fold response cross-validation. All category
indicators of a held-out subject-variable outcome are excluded together.
The same balanced folds (seed 20261019) are used at every penalty, preserving
the caller's random state. Fits are warm-started along decreasing penalties.
Each candidate needs five fits; the winner is then refitted on all observations.
The criterion is average held-out negative log probability. A nonconverged
fit stops the run; increase the iteration limit before retrying.
Ties favor the larger penalty.
The selected penalty is fixed in all null and bootstrap fits. Rank-zero observed
fits return empty signal spaces. Missing bootstrap directions receive a
90-degree perturbation bound and zero projection-cache columns. Zero-rank null
fits have undefined (NA) angles; null-angle culling is unavailable if no valid
calibration remains. Diagnostics record the actual refit ranks and valid null draws.

The null model sets theta equal to mu for each subject. Fit one independent
null reference subspace and compare each of null_nsim further draws' leading
directions to that fixed subspace using the orthogonal projection length.
For MCA this retains the previous fixed-rank SVD procedure. For nuclear norm,
only nonzero fitted directions are used, up to the observed rank.
}
