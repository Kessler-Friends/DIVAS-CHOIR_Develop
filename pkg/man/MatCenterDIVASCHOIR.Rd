\name{MatCenterDIVASCHOIR}
\alias{MatCenterDIVASCHOIR}
\title{Construct MCA input for categorical DIVAS blocks}
\usage{
MatCenterDIVASCHOIR(block, allow_empty_levels = FALSE, method = c("mca", "nuclear_norm"))
}
\arguments{
\item{block}{A data frame or matrix with subjects in rows and categorical variables in columns.}
\item{method}{Preparation method. MCA centers and scales one-hot indicators; nuclear_norm preserves zero/one observations and enforces latent centering during likelihood estimation.}
\item{allow_empty_levels}{If TRUE, retain an unobserved factor level as a zero row. Bootstrap draws use this to keep the same category rows as the observed block.}
}
\value{
A method-tagged category-level-by-subject matrix with marginal-probability and level metadata.
}
\description{
Encodes observations and applies the selected preparation method. The CHOIR extractor owns the method-specific preparation for original and bootstrap data.
}
