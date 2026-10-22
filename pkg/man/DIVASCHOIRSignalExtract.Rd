\name{DIVASCHOIRSignalExtract}
\alias{DIVASCHOIRSignalExtract}
\title{Categorical signal extraction for DIVAS-CHOIR}
\usage{
DIVASCHOIRSignalExtract(datablock, nsim = 400, iplot = FALSE,
  cull = 0.382, seed = NULL, fixed_rank = NULL, lambda = NULL,
  cull_rank = FALSE)
}
\arguments{
\item{datablock}{A list of centered and scaled MCA matrices, with category levels in rows and subjects in columns. Use \code{MatCenterDIVASCHOIR} to prepare each block.}
\item{nsim}{Number of bootstrap samples.}
\item{iplot}{Whether to plot singular values.}
\item{cull}{Multiplier applied to the fifth percentile of the simulated categorical null angles for rank culling.}
\item{seed}{Optional random seed.}
\item{fixed_rank}{Optional initial rank, either shared or one value per block, before culling.}
\item{lambda}{Optional positive nuclear-norm penalty, either shared or one value per block. Singular values are soft-thresholded by this value.}
\item{cull_rank}{Apply rank culling using simulated categorical null angles. Defaults to FALSE, retaining all initial signal directions.}
}
\value{
A list containing categorical signal bases, ranks, singular values, bootstrap angle bounds, and \code{nullAngleCaches}. Each null cache contains 1,000 sample-space (\code{subject}) and category-space (\code{feature}) angles in degrees.
}
\description{
Computes categorical MCA signal spaces and parametric-bootstrap angle bounds.
}
\details{
The null model sets theta equal to the category intercept mu for every sample.
Softmax of mu gives the observed category proportions within each variable.
Simulate one reference categorical dataset, recompute its MCA matrix, and
estimate its rank-r SVD subspaces at the initial fitted rank. Hold these
subspaces fixed across 1,000 further independent null datasets. The angle from
each new draw's leading estimated direction to the fixed subspace defines the sample-space and category-space
null distributions. It is the arccosine of the length of the unit direction's
orthogonal projection onto that subspace, in degrees. Their fifth
percentiles replace the uniform random-direction references.
}
