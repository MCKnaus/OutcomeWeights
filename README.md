# OutcomeWeights   <img src="man/figures/logo.png" align="right" alt="OutcomeWeights logo" height="150"/>

<!-- badges: start -->
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/OutcomeWeights)](https://CRAN.R-project.org/package=OutcomeWeights)
[![License](https://img.shields.io/badge/license-GPL--3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html)
[![Downloads_Total](https://cranlogs.r-pkg.org/badges/grand-total/OutcomeWeights)](https://CRAN.R-project.org/package=OutcomeWeights)
[![Downloads_Monthly](https://cranlogs.r-pkg.org/badges/OutcomeWeights)](https://CRAN.R-project.org/package=OutcomeWeights)
[![Project_Status](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
<!-- badges: end -->

This R package calculates the outcome weights of [Knaus (2024)](https://arxiv.org/abs/2411.11559). 
Its use is illustrated in the vignettes of the [package website](https://mcknaus.github.io/OutcomeWeights/).

The core functionality is the `get_outcome_weights()` method implementing the theoretical result in Proposition 1 of the paper. It shows that the outcome weights vector can be obtained in the general form
$\boldsymbol{\omega'} = (\boldsymbol{\tilde{Z}'\tilde{D}})^{-1} \boldsymbol{\tilde{Z}'T}$
where $\boldsymbol{\tilde{Z}}$, $\boldsymbol{\tilde{D}}$ and $\boldsymbol{T}$ are pseudo-instrument, pseudo-treatment and the transformation matrix, respectively. 

In the future it should be compatible with as many estimated R objects as possible.

The package can be downloaded from CRAN:
```R
install.packages("OutcomeWeights")
```

The package is work in progress. Find here the current state (suggestions welcome):

### In progress
- [ ] Compatibility with [`grf`](https://grf-labs.github.io/grf/) package
  - [x] `causal_forest()` outcome weights for CATE
  - [x] `instrumental_forest()` outcome weights CLATE
  - [x] `causal_forest()` outcome weights for ATE from `average_treatment_effect()`
  - [ ] All outcome weights for average parameters compatible with `average_treatment_effect()`
- [ ] Package internal Double ML implementation handling the required outcome smoother matrices
  - [x] Nuisance parameter estimation based on honest random forest (`regression_forest()` of `grf` package)
  - [x] `dml_with_smoother()` function runs for PLR, PLR-IV, AIPW-ATE, and Wald_AIPW and is compatible with `get_outcome_weights()`
  - [ ] Add more Double ML estimators
  - [ ] Add support for more smoothers

### Envisioned features
- [ ] Compatibility with [`DoubleML`](https://docs.doubleml.org/stable/index.html) (this is a non-trivial task as the `mlr3` environment it builds on does not provide smoother matrices)
  - [ ] Extract the smoother matrices of `mlr3` available, where possible
  - [ ] Make the smoother matrices of `mlr3` accessible within DoubleML
  - [ ] Write `get_outcome_weights()` method for DoubleML estimators
- [ ] Collect packages where weights could be extracted and implement them


The development version is available using the `devtools` package:
```R
library(devtools)
install_github(repo="MCKnaus/OutcomeWeights")
```

## References

Knaus, M. C. (2024). Treatment effect estimators as weighted outcomes, [arXiv:2411.11559](https://arxiv.org/abs/2411.11559)
