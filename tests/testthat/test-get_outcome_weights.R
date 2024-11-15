test_that("get_outcome_weights computes correct weights", {
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
  
  # Run DML and get weights
  dml = dml_with_smoother(Y,D,X,Z,n_cf_folds = 2)
  results_dml = summary(dml)
  omega_dml = get_outcome_weights(dml)
  summary(omega_dml)
  
  # Run causal/instrumental forest with a pre-fitted outcome smoother
  forest.Y = grf::regression_forest(X, Y)
  Y.hat = predict(forest.Y)$predictions
  outcome_smoother = grf::get_forest_weights(forest.Y)
  c.forest = grf::causal_forest(X, Y, D, Y.hat = Y.hat)
  cates = predict(c.forest)$predictions
  i.forest = grf::instrumental_forest(X, Y, D, Z, Y.hat = Y.hat)
  clates = predict(i.forest)$predictions
  
  # Calculate outcome weights
  omega_cf = get_outcome_weights(c.forest, S = outcome_smoother)
  omega_if = get_outcome_weights(i.forest, S = outcome_smoother)
  
  # Expectation about weights replicating point estimates
  expect_equal(as.numeric(omega_dml$omega %*% Y), 
               as.numeric(results_dml[,1]))
  expect_equal(as.numeric(omega_cf$omega %*% Y), 
               as.numeric(cates))
  expect_equal(as.numeric(omega_if$omega %*% Y), 
               as.numeric(clates))
})
