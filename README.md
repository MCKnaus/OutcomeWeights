# Outcome Weights

<!-- badges: start -->
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/OutcomeWeights)](https://CRAN.R-project.org/package=OutcomeWeights)
[![License](https://img.shields.io/badge/license-GPL--3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html)
[![Downloads_Total](https://cranlogs.r-pkg.org/badges/grand-total/OutcomeWeights)](https://CRAN.R-project.org/package=OutcomeWeights)
[![Downloads_Monthly](https://cranlogs.r-pkg.org/badges/OutcomeWeights)](https://CRAN.R-project.org/package=OutcomeWeights)
[![Project_Status](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
<!-- badges: end -->

This R package calculates the outcome weights of [Knaus (2024)](https://arxiv.org/abs/2411.11559). Its use is illustrated in the 
[average effects R notebook](https://mcknaus.github.io/assets/code/Notebook_Application_average_401k.nb.html) and the
[heterogeneous effects R notebook](https://mcknaus.github.io/assets/code/Notebook_Application_heterogeneous_401k.nb.html) 
as supplementary material to the paper.

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


The following code creates synthetic data to showcase how causal forest weights are extracted and that they perfectly replicate the original output:

``` r
if (!require("OutcomeWeights")) install.packages("OutcomeWeights", dependencies = TRUE)
library(OutcomeWeights)

# Sample from DGP borrowed from grf documentation
n = 500
p = 10
X = matrix(rnorm(n * p), n, p)
W = rbinom(n, 1, 0.5)
Y = pmax(X[, 1], 0) * W + X[, 2] + pmin(X[, 3], 0) + rnorm(n)

# Run outcome regression and extract smoother matrix
forest.Y = grf::regression_forest(X, Y)
Y.hat = predict(forest.Y)$predictions
outcome_smoother = grf::get_forest_weights(forest.Y)

# Run causal forest with external Y.hats
c.forest = grf::causal_forest(X, Y, W, Y.hat = Y.hat)

# Predict on out-of-bag training samples.
cate.oob = predict(c.forest)$predictions

# Predict using the forest.
X.test = matrix(0, 101, p)
X.test[, 1] = seq(-2, 2, length.out = 101)
cate.test = predict(c.forest, X.test)$predictions

# Calculate outcome weights
omega_oob = get_outcome_weights(c.forest, S = outcome_smoother)
omega_test = get_outcome_weights(c.forest, S = outcome_smoother, newdata = X.test)

# Observe that they perfectly replicate the original CATEs
all.equal(as.numeric(omega_oob$omega %*% Y), 
          as.numeric(cate.oob))
all.equal(as.numeric(omega_test$omega %*% Y), 
          as.numeric(cate.test))

# Also the ATE estimates are perfectly replicated
omega_ate = get_outcome_weights(c.forest,target = "ATE", S = outcome_smoother, S.tau = omega_oob$omega)
all.equal(as.numeric(omega_ate$omega %*% Y),
          as.numeric(grf::average_treatment_effect(c.forest, target.sample = "all")[1]))
```


The development version is available using the `devtools` package:
```R
library(devtools)
install_github(repo="MCKnaus/OutcomeWeights")
```

## References

Knaus, M. C. (2024). Treatment effect estimators as weighted outcomes, [arXiv:2411.11559](https://arxiv.org/abs/2411.11559)
