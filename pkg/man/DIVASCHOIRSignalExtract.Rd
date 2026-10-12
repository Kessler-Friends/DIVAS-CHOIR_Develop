\name{DIVASCHOIRSignalExtract}
\alias{DIVASCHOIRSignalExtract}
\title{Categorical signal extraction for DIVAS-CHOIR}
\usage{
DIVASCHOIRSignalExtract(datablock, nsim = 400, iplot = FALSE,
  cull = 0.382, seed = NULL, fixed_rank = NULL, lambda = NULL,
  cull_rank = TRUE)
}
\arguments{
\item{datablock}{A list of centered and scaled MCA matrices, with category levels in rows and subjects in columns. Use \code{MatCenterDIVASCHOIR} to prepare each block.}
\item{nsim}{Number of bootstrap samples.}
\item{iplot}{Whether to plot singular values.}
\item{cull}{Multiplier applied to the random-angle reference for rank culling.}
\item{seed}{Optional random seed.}
\item{fixed_rank}{Optional initial rank, either shared or one value per block, before culling.}
\item{lambda}{Optional positive nuclear-norm penalty, either shared or one value per block. Singular values are soft-thresholded by this value.}
\item{cull_rank}{Apply the original DIVAS random-angle rank culling rule.}
}
\value{
A list containing categorical signal bases, ranks, singular values, and bootstrap angle bounds.
}
\description{
Computes categorical MCA signal spaces and parametric-bootstrap angle bounds.
}
