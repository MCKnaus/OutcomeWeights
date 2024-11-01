###### The get_outcome_weights method ###############
## this part is needed for using an get_outcome_weights.whatever method
## methods for different packages are provided in separate <package_name>.R files.

#' Outcome weights method
#'
#' @description This is a generic method for getting outcome weights. 
#' It calculates the outcome weights for objects created by other packages.
#' See get_outcome_weight.<compatible_fct> in the package documentation for compatible functions.
#'
#' @param object An object, obtained from other packages.
#' @param ... Additional arguments specific to object class implementations. 
#' See the documentation which object requires which additional arguments.
#'
#' @return A vector or matrix of outcome weights.
#'
#' @export
#'
get_outcome_weights = function(object,...) UseMethod("get_outcome_weights") # makes get_outcome_weights a generic method


#' Outcome weights maker for pseudo-IV estimators.
#'
#' @description This is a generic function taking pseudo-instrument,
#' pseudo-treatment and the transfmoration matrix as inputs and returning 
#' outcome weights
#'
#' @param Z.tilde Numeric vector of pseudo-instrument outcomes.
#' @param D.tilde Numeric vector of pseudo-treatment.
#' @param T_mat Transformation matrix
#'
#' @return A vector of outcome weights.
#'
#' @export
#'
pive_weight_maker = function(Z.tilde,D.tilde,T_mat) {
  omega = (t(Z.tilde) %*% D.tilde)^(-1) %*% t(Z.tilde) %*% T_mat
  return(omega)
}

