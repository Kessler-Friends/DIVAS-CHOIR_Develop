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

DIVAS requires `CVXR (>= 1.8.2)`. Nuclear multinomial fits use its declarative
objective and constraints with the Clarabel solver.

```R
# Install devtools (if not already installed)
install.packages("devtools")

# Install CVXR 1.8.2 or newer.
install.packages("CVXR")
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

### Choosing the CHOIR estimator

`DIVASmain()` accepts raw categorical data frames and passes
`categorical_method = "mca"` (the default) or `"nuclear_norm"` to the CHOIR
extractor. Continuous estimation is unchanged.

`pkg/R/SignalExtractUtils.R` contains shared preparation in
`MatCenterDIVASCHOIR()`, MCA preparation in `choir_prepare_mca()`, and nuclear
setup and latent centering in `choir_prepare_nuclear()`.
`pkg/R/DIVASCHOIREstimation.R` contains the MCA and nuclear estimators.
`pkg/R/DIVASCHOIRSignalExtract.R` selects the method and handles observed fits,
null simulations, bootstrap fits, angles and outputs.

- MCA: one-hot encode, subtract empirical category proportions, and scale by
  their square roots; estimate the signal with truncated SVD.
- Nuclear: leave one-hot observations unchanged and initialize category groups,
  the centered intercept and centered latent eta; fit multinomial likelihood
  plus a nuclear-norm penalty on the latent signal. Center eta across subjects
  and enforce within-variable sum-to-zero constraints during estimation.

The extractor calls the same method for the observed fit, null reference fits,
and both bootstrap stages, including rebuilding the correct input each time.
Nuclear probabilities come directly from fitted intercept-plus-eta logits;
the fit is expressed as multinomial negative log likelihood plus
`lambda * CVXR::norm_nuc(eta)`, with linear centering constraints.
there is no MCA back-map in that path. Downstream reconstruction uses the MCA
matrix for MCA or the fitted latent eta for nuclear. Detailed results retain
both prepared observations and the analysis matrices.

Nuclear lambda selection uses five-fold response CV. Each held-out unit is one
subject-variable outcome, with all its category indicators excluded together.
Choose the lowest average held-out log loss, then fit that penalty on all data.
`categorical_rank` applies only to MCA. Folds use a fixed seed (20261019), preserve
the caller's random state, and retain training observations for each subject and
variable. Warm starts follow decreasing penalties within each fold.

`categorical_lambda_grid` supplies candidate nonnegative penalties. The default
uses factors 1.05, 0.75, 0.45, 0.25, and 0.1 times the intercept-only score spectral
norm. Thus the default searches a finite positive grid, including an intercept-only
candidate, rather than optimizing over all nonnegative lambda. Selection needs
five fits per candidate and one final full-data fit. A nonconverged fit stops
the run so the iteration limit can be increased; ties favor
the larger lambda. The score table is saved in `Phase1_categorical$estimation`.

The chosen lambda stays fixed in null and bootstrap refits. Their ranks may vary.
Absent bootstrap directions receive 90-degree bounds and zero projection-cache
columns. Null fits with rank zero have NA angles, not arbitrary directions;
diagnostics record actual ranks and valid draws. Optional null-angle culling
stops if no valid calibration remains. An observed rank-zero fit returns empty
spaces without running null simulations or bootstraps.

```R
result <- DIVASmain(
  continuous_datablock = continuous_blocks,
  categorical_datablock = categorical_blocks,
  categorical_method = "nuclear_norm",
  categorical_lambda_grid = c(0.5, 1, 2, 4), # omit to use the default grid
  categorical_cull = FALSE,
  ReturnDetail = TRUE
)
```

`Phase1_categorical$estimation` records method, rank, penalty and bootstrap
metadata. Existing MCA preparation calls still work; a matrix prepared for one
method cannot be reused as input to the other method. Pass raw data instead.

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


### Categorical reconstruction

After joint directions are selected, `DJIVEReconstructMJ()` routes continuous
blocks to the existing least-squares `MatReconstructMJ()` and categorical blocks
to `choir_reconstruct()` in `R/DIVASCHOIRReconstruct.R`. This applies to both MCA
and nuclear extraction. Reconstruction uses the original one-hot responses,
not the standardized MCA matrix or a least-squares projection of estimated eta.

For each categorical block, all retained joint, partially shared and individual
score directions are fixed and fitted together. The multinomial likelihood
estimates an unpenalized intercept and category coefficients. Scores are centered
across subjects; category sums of intercepts and coefficients are constrained to
zero for each variable. Components add in logit space; probabilities follow from
softmax within each variable after adding the intercept. No new penalty or rank
selection is performed. Solver failures stop explicitly; extreme logits warn of
possible separation and the lack of a finite unpenalized optimum.

`matLoadings` and categorical `Loadings` contain fitted category coefficients;
`matBlocks` contain their additive latent effects. Categorical `scoresList` keeps
the centered fixed score coordinates, without the continuous inverse-loading
transformation. `categoricalReconstruction` (also included in compact results)
provides `mu`, `eta`, `theta`, `probabilities`, `residual` (one-hot minus fitted
probabilities), the centered `basis`, and likelihood/solver diagnostics. Historical
saved fits keep their previous reconstruction until explicitly refitted.
