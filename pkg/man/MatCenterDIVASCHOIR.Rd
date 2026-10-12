\name{MatCenterDIVASCHOIR}
\alias{MatCenterDIVASCHOIR}
\title{Construct MCA input for categorical DIVAS blocks}
\usage{
MatCenterDIVASCHOIR(block, allow_empty_levels = FALSE)
}
\arguments{
\item{block}{A data frame or matrix with subjects in rows and categorical variables in columns.}
\item{allow_empty_levels}{If TRUE, retain an unobserved factor level as a zero row. Bootstrap draws use this to keep the same category rows as the observed block.}
}
\value{
A centered and scaled category-level-by-subject matrix with marginal-probability and level metadata.
}
\description{
One-hot encodes the categorical variables and forms the weighted MCA SVD input.
}
