# DIVAS: Data Integration via Analysis of Subspaces

<!-- badges: start -->
[![R-CMD-check](https://github.com/ByronSyun/DIVAS_Develop/workflows/R-CMD-check/badge.svg)](https://github.com/ByronSyun/DIVAS_Develop/actions)
<!-- badges: end -->

<p align="center">
<img src="man/figures/DIVAS_logo.png" width="200" alt="DIVAS Logo">
</p>

## Introduction

We provide an R implementation of DIVAS (Prothero et al., 2024), a statistical method for multi-modal data integration. Via statistical analysis of subspaces, DIVAS identifies joint, partially shared, and individual variation across multiple data blocks in a completely data-driven way. In addition to translating the original MATLAB implementation into an accessible R package, we also provide visualization tools and worked examples for exploring DIVAS results in practice.

**Documentation website**: [https://byronsyun.github.io/DIVAS_Develop/](https://byronsyun.github.io/DIVAS_Develop/)

## Installation

### Dependencies

The DIVAS package requires the current 1.x line of the `CVXR` package for compatibility with the SCS solver interface. The package has been tested with `CVXR` 1.0-15. In `DESCRIPTION`, DIVAS also declares `CVXR (>= 1.0-15)` as a package dependency.

```R
# Install devtools (if not already installed)
install.packages("devtools")

# Install CVXR 1.x. DIVAS has been tested with CVXR 1.0-15.
install.packages("remotes")
remotes::install_version("CVXR", version = "1.0-15", repos = "https://cloud.r-project.org")

# Alternatively, install the latest CRAN release if it is compatible with your R version.
# install.packages("CVXR")
```

### Installing the DIVAS package

You can install the development version of DIVAS from GitHub using `devtools`:

```R
# Install DIVAS package from the main branch on GitHub
devtools::install_github("ByronSyun/DIVAS_Develop/pkg", ref = "main")

# Or install from a local folder if you have cloned the repository
# devtools::install("path/to/DIVAS-main/pkg")
```

## Usage Examples

The DIVAS package supports analysis of various data formats. Here is a simple example using the built-in toy MATLAB dataset:

```R
library(R.matlab)
library(DIVAS)

data_path <- system.file("extdata", "toyDataThreeWay.mat", package = "DIVAS")
data <- readMat(data_path)

datablock <- list(
  X1 = data$datablock[1,1][[1]][[1]],
  X2 = data$datablock[1,2][[1]][[1]],
  X3 = data$datablock[1,3][[1]][[1]]
)

result <- DIVASmain(continuous_datablock = datablock)
dataname <- paste0("DataBlock_", 1:length(datablock))
plots <- DJIVEAngleDiagnosticJP(datablock, dataname, result, 566, "Demo")
print(plots)
```

For more detailed tutorials, see the documentation website and linked case studies below.

### Categorical nuclear-norm rank selection

The optional `nuclear_cv` method selects a penalty separately for each
categorical block. It holds out complete subject-by-variable responses, fits a
nuclear-norm-penalized MCA matrix on the remaining responses, and chooses the
penalty with the lowest held-out squared error in MCA coordinates. The full
categorical block is then refit by soft-thresholding its singular values.
The default retains the initial categorical rank (set with
`categorical_rank`, defaulting to at most five) and uses multinomial
parametric bootstrap samples to estimate angle bounds. Set
`categorical_cull = TRUE` to apply angle culling explicitly. The reference
distribution uses categorical datasets simulated under
`theta = mu`: softmax of the category intercepts gives the observed proportions,
One null draw is centered/scaled as MCA and supplies fixed rank-r reference
subspaces. Each of 1,000 further independent null draws supplies a leading
direction. The fifth percentile of its angles to the fixed subspace sets the
reference in both sample and category spaces. The simulated angles
are in `result$Phase1_categorical$nullAngleCaches`.
The retained rank and cutoffs are in
`result$Phase1_categorical$rankDiagnostics`.
Bootstrap MCA uses empirical category proportions without an added count.
If a level is absent in a simulated draw, its row is zero in the common
category space. Nuclear cross-validation likewise omits levels absent from
a training fold when computing that fold's validation error.

```R
result <- DIVASmain(
  continuous_datablock = list(RNA = rna_matrix),
  categorical_datablock = list(clinical = clinical_data_frame),
  categorical_method = "nuclear_cv",
  categorical_cv_folds = 5,
  seed = 42
)
result$categorical_lambda
result$categorical_cv[[1]]$cv
```

`clinical_data_frame` has one row per subject and one column per categorical
variable. `rna_matrix` has features in rows and the same subjects in columns.
The selected penalty is in the scale of the centered, weighted MCA matrix;
the CV score is a reconstruction error, rather than multinomial likelihood.
Use `categorical_lambda_grid` to supply candidate penalties explicitly.

## Available Datasets

We provide the following examples to illustrate the use of DIVAS in different scenarios.

| Dataset             | Brief Description                                  | Vignette Link                                                                                              | Format | Primary Reference      |
|---------------------|----------------------------------------------------|------------------------------------------------------------------------------------------------------------|--------|------------------------|
| toyDataThreeWay.mat | Synthetic 3-block data with known joint structures | [Toy Dataset Example](https://byronsyun.github.io/DIVAS_Develop/articles/DIVAS_Toy_Dataset_Example.html) | .mat   | Prothero et al. (2024) |
| toyDataThreeWayCategorical.rds | Mixed version of the three-block toy: X1 and X2 continuous, X3 converted to 36 binary traits from its three leading sample directions | `data-raw/build_toy_threeway_categorical.R` | .rds | Derived from Prothero et al. (2024) |
| gnp_imputed.qs      | GNP economic time series data                      | [GNP Dataset Example](https://byronsyun.github.io/DIVAS_Develop/articles/DIVAS_GNP_Dataset_Example.html) | .qs    | Stock & Watson (2016)  |
| COVID-19 Multi-Omics | 6-block integration: scRNA-seq (4 cell types), proteomics, metabolomics from 114 COVID-19 patient samples | [COVID Case Study](https://byronsyun.github.io/DIVAS_COVID19_CaseStudy/) | .rds   | Su et al. (2020) |

The mixed toy can be loaded with
`readRDS(system.file("extdata", "toyDataThreeWayCategorical.rds", package = "DIVAS"))`.
Pass its `continuous` and `categorical` lists to `DIVASmain()`.
The conversion is deterministic and uses X3 alone; the R script in `data-raw`
regenerates it from the original MATLAB file.

## Case Study: COVID-19 Multi-Omics Analysis

This project serves as a comprehensive, real-world application of the DIVAS package on a complex multi-omics dataset from a COVID-19 patient cohort. It demonstrates the full data processing and analysis workflow, from raw data cleaning to final DIVAS results, showcasing the practical utility of the package.

**➡️ [View the full analysis on GitHub](https://github.com/ByronSyun/DIVAS_COVID19_CaseStudy)**

## Developers

* **[Jiadong Mao](https://github.com/jiadongm)** - *Lead Developer, Maintainer*
* **[Yinuo Sun](https://github.com/ByronSyun)** - *Package Developer, Maintainer*

## References

Prothero, J., et al. (2024). Data integration via analysis of subspaces (DIVAS).

Su, Y., Chen, D., Yuan, D., et al. (2020). Multi-Omics Resolves a Sharp Disease-State Shift between Mild and Moderate COVID-19. Cell, 183(6), 1479-1495. https://doi.org/10.1016/j.cell.2020.10.037

## License

This project is licensed under the GNU Affero General Public License v3.0 (AGPL-3) - see the LICENSE file for details.
