#' Calls C++ implementation to calculate standardized mean differences.
#' 
#' Calculates standardized mean differences between treated and controls and towards target means 
#' for an outcome weights matrix with potentially many rows like for CATEs.
#' 
#' @param X Covariate matrix with N rows and p columns.
#' @param treat Binary treatment variable.
#' @param omega Outcome weights matrix with dimension number of weight vectors for which balancing should be checked 
#' x number of training units.
#' @param target Optional matrix with dimension number of weight vectors for which balancing should be checked
#' x p indicating the target values the covariates should be balanced towards. 
#' If NULL, average of X used as target of ATE.
#' 
#' @return 3D-array of dimension p x 6 x number of weight vectors for which balancing should be checked 
#' where the second dimension provides the following quantities:
#' - "Mean 0": The weighted control mean
#' - "Mean 1": The weighted treated mean
#' - "SMD balancing": Standardized mean differences for covariate balancing (Mean 1 - Mean 0) / sd(X) 
#' - "SMD targeting 0": Standardized mean difference to assess targeting of control (Mean 0 - target) / sd(X) 
#' - "SMD targeting 1": Standardized mean difference to assess targeting of treated (Mean 1 - target) / sd(X) 
#' 
#' @export
#' 
standardized_mean_differences = function(X,treat,omega,
                                         target=NULL) {

  if (!is.matrix(X) | !is.matrix(omega) | !is.matrix(target)) stop("X, omega and target must be matrices.")
  if (nrow(omega) != nrow(target)) stop("Please provide omega and target matrix with same number of rows.")
  if (is.null(colnames(X))) stop("Please provide covariate matrix X with column names.")
  if (is.null(target)) target = matrix(rep(colMeans(X), nrow(omega)), nrow = nrow(omega), byrow = TRUE)
  if (ncol(X) != ncol(target)) stop("Please provide covariate matrix X with column names.")
  
  output = smd_rcpp(X,treat,omega,target)
  dimnames(output) = list(colnames(X),
                          c("Mean 0","Mean 1","SMD balancing","SMD targeting 0","SMD targeting 1"),
                          seq(1:nrow(omega)))
  
  class(output) = c("standardized_mean_differences")
  
  return(output)
}



#' \code{summary} method for class \code{\link{standardized_mean_differences}}
#'
#' Calls a C++ function to quickly summarize potentially many standardized mean differences.
#' 
#' @param object Object of class \code{\link{standardized_mean_differences}}.
#' @param ... further arguments passed to summary method.
#' 
#' @return 3D-array of dimension 
#' - c("Maximum absolute SMD","Mean absolute SMD","Median absolute SMD", / % of absolute SMD > 20",
#' "# / % of absolute SMD > 10","# / % of absolute SMD > 5") x
#' - c("Balancing","Targeting") x 
#' - number of weight vectors for which balancing should be checked 
#' 
#' @export
#' 
summary.standardized_mean_differences = function(object, ...) {
  
  output = summary_smd_rcpp(object)
  dimnames(output) = list(c("Maximum absolute SMD","Mean absolute SMD","Median absolute SMD",
                            "# / % of absolute SMD > 20","# / % of absolute SMD > 10","# / % of absolute SMD > 5"),
                          c("Balancing","Targeting"),
                          seq(1:dim(output)[3]))
  
  return(output)
}

