test_that("summary_weights has correct dimensions", {
  # Set up sample data
  set.seed(1234)
  n = 100
  p = 3
  X = matrix(rbinom(n * p, 1, 0.5), n, p)
  Z = rbinom(n, 1, 0.5)
  Q = rbinom(n, 1, 0.5)
  D = Q * Z
  tau = X[, 1] / 2
  Y = rowSums(X[, 1:3]) + tau * D + Q + rnorm(n)
  
  # Run causal/instrumental forest with a pre-fitted outcome smoother
  forest.Y = grf::regression_forest(X, Y)
  Y.hat = predict(forest.Y)$predictions
  outcome_smoother = grf::get_forest_weights(forest.Y)
  c.forest = grf::causal_forest(X, Y, D, Y.hat = Y.hat)
  
  # Calculate outcome weights
  omega_cf = get_outcome_weights(c.forest, S = outcome_smoother)
  
  # Get summary weights
  summary_weights = summary(omega_cf,quiet = T)
  
  # Check SMDs
  cov_bal = standardized_mean_differences(X,D,omega_cf$omega,X)
  sum_cov_bal = summary(cov_bal)
  
  # Check if summary_weights has correct dimensions
  expect_equal(dim(summary_weights), c(2, n, 6))
  expect_equal(dim(cov_bal), c(p, 5, n))
  expect_equal(dim(sum_cov_bal), c(6, 2, n))
})