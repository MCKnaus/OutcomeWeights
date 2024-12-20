test_that("all estimator combinations run", {
  skip_on_cran()
  # Set up sample data
  set.seed(1234)
  n = 100
  p = 3
  X = matrix(rbinom(n * p, 1, 0.5), n, p)
  Z = rbinom(n, 1, 0.5)
  Q = rbinom(n, 1, 0.5)
  D = Q * Z
  tau =  X[, 1] / 2
  Y = rowSums(X[, 1:3]) + tau * D + Q + rnorm(n)
  
  # Run all possible combinations with the minimum required variables
  dml1 = dml_with_smoother(Y,D,X,n_cf_folds = 2,
                          estimators = "PLR")
  dml2 = dml_with_smoother(Y,D,X,Z,n_cf_folds = 2,
                           estimators = "PLR_IV")
  dml3 = dml_with_smoother(Y,D,X,n_cf_folds = 2,
                           estimators = "AIPW_ATE")
  dml4 = dml_with_smoother(Y,D,X,Z,n_cf_folds = 2,
                           estimators = "Wald_AIPW")
  dml5 = dml_with_smoother(Y,D,X,Z,n_cf_folds = 2,
                           estimators = c("PLR","PLR_IV"))
  dml6 = dml_with_smoother(Y,D,X,n_cf_folds = 2,
                           estimators = c("PLR","AIPW_ATE"))
  dml7 = dml_with_smoother(Y,D,X,Z,n_cf_folds = 2,
                           estimators = c("PLR","Wald_AIPW"))
  dml8 = dml_with_smoother(Y,D,X,Z,n_cf_folds = 2,
                           estimators = c("AIPW_ATE","Wald_AIPW"))
  dml9 = dml_with_smoother(Y,D,X,Z,n_cf_folds = 2,
                           estimators = c("PLR", "PLR_IV", "AIPW_ATE"))
  dml10 = dml_with_smoother(Y,D,X,Z,n_cf_folds = 2,
                           estimators = c("PLR", "PLR_IV", "Wald_AIPW"))
  dml11 = dml_with_smoother(Y,D,X,Z,n_cf_folds = 2,
                           estimators = c("PLR", "AIPW_ATE", "Wald_AIPW"))

  expect_equal(class(dml1), "dml_with_smoother")
  expect_equal(class(dml2), "dml_with_smoother")
  expect_equal(class(dml3), "dml_with_smoother")
  expect_equal(class(dml4), "dml_with_smoother")
  expect_equal(class(dml5), "dml_with_smoother")
  expect_equal(class(dml6), "dml_with_smoother")
  expect_equal(class(dml7), "dml_with_smoother")
  expect_equal(class(dml8), "dml_with_smoother")
  expect_equal(class(dml9), "dml_with_smoother")
  expect_equal(class(dml10), "dml_with_smoother")
  expect_equal(class(dml11), "dml_with_smoother")
})
