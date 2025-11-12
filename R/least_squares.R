#' Outcome weights for the \code{\link[estimatr]{lm_robust}} command
#'
#' @description Post-estimation command to calculate the outcome weights for linear regression implemented via the \code{\link[estimatr]{lm_robust}} command.
#'
#' @param object An object of class \code{lm_robust}, i.e. the result of running \code{\link[estimatr]{lm_robust}}.
#'
#' @param treat_pos The position of the treatment indicator in the model matrix. By default 2,
#' which is appropriate for the case with an intercept and the treatment indicator as the first explanatory variable.
#'
#' @param ... Pass potentially generic \link{get_outcome_weights} options.
#'
#' @return \link{get_outcome_weights} object with `omega` containing weights and `treat` the treatment
#'
#' @export
#'
get_outcome_weights.lm_robust = function(object,treat_pos=2,...) {
  n = object$nobs
  weighted = object$weighted
  if (weighted){
    ncol = dim(as.matrix(model.frame(object)))[2]
    temp = as.matrix(model.frame(object))[,c(2:(ncol-1))]# the weights are in the last column, Y is in the first
  } else{
    temp = as.matrix(model.frame(object))[,-1] # no need to delete a column for the weights
  }
  
  # now build the matrix of regressors from the model frame.
  # This entails checking whether there was an intercept in the original call
  if(dim(temp)[2]==length(object$coefficients)){ # if this condition is true, no intercept
    x = as.matrix(temp[,!is.na(object$coefficients)])
  } else if (dim(temp)[2]+1==length(object$coefficients)){ # if this condition is true, add intercept
    temp = cbind(rep(1,n),temp) # add the intercept
    x = as.matrix(temp[,!is.na(object$coefficients)]) # use only the columns with coefficients
  } else {warning("something went wrong")}
  
  # if there are NA coefficients, we may need to extract a different row of the matrix
  # (as we deleted the cols corresponding to the NA coefficients)
  
  C = matrix(0,1,ncol(x))
  
  if (sum(is.na(object$coefficients))==0){ # everything is fine
    C[treat_pos] = 1
  } else if (sum(which(is.na(object$coefficients))<=treat_pos)==0){ # if the NA coefficients are after the TE -> no problem
    C[treat_pos] = 1
  } else if (sum(which(is.na(object$coefficients))<=treat_pos)>0){
    num_NA = sum(which(is.na(object$coefficients))<=treat_pos) # number of NA coefficients -> del. cols before TE
    C[treat_pos-num_NA] = 1 # if some cols were deleted before the TE, we need a different row of the resp. matrix
  } else {stop("Check for multicollinearity")} # just in case something is wrong -> but shouldn't be the case
  
  D = x[,treat_pos]
  
  ### Calculation of the outcome weights
  # if sample weights were included in the original call, weights should take care of this
  if (weighted){
    weights = object$weights
    x = x*sqrt(weights)
    omega = C %*% solve(crossprod(x),tol=.Machine$double.xmin) %*% t(x)*sqrt(weights)
  } else{
    omega = C %*% solve(crossprod(x),tol=.Machine$double.xmin) %*% t(x)
  }
  # check whether OLS coefficient and omega*y are equal, warning if not
  Y = as.matrix(model.frame(object)[,1],n)
  if(!isTRUE(all.equal(as.numeric(omega %*% Y),
                       as.numeric(object$coefficients[treat_pos])))){
    warning("Estimated Treatment Effects using weights differ from original estimates.")
  }
  
  output = list(
    "omega" = matrix(omega,nrow=1),
    "treat" = D
  )
  
  class(output) = c("get_outcome_weights")
  
  return(output)
}


#' Outcome weights for the lm command
#'
#' @description Post-estimation command to calculate the outcome weights for linear regression implemented via the \code{\link[stats]{lm}} command.
#'
#' @param object An object of class \code{lm}, i.e. the result of running \code{\link[stats]{lm}}.
#'
#' @param treat_pos The position of the treatment indicator in the model matrix. By default 2,
#' which is appropriate for the case with an intercept and the treatment indicator as the first explanatory variable.
#'
#' @param ... Pass potentially generic \link{get_outcome_weights} options.
#'
#' @return \link{get_outcome_weights} object with `omega` containing weights and `treat` the treatment
#'
#' @export
#'
get_outcome_weights.lm = function(object,treat_pos=2,...) {
  n = length(object$residuals)
  temp = model.matrix(object)
  x = temp[,!is.na(object$coefficients)]
  C = matrix(0,1,ncol(x))
  D = x[,treat_pos]
  
  if (sum(is.na(object$coefficients))==0){ # everything is fine
    C[treat_pos] = 1
  } else if (sum(which(is.na(object$coefficients))<=treat_pos)==0){ # if the NA coefficients are after the TE -> no problem
    C[treat_pos] = 1
  } else if (sum(which(is.na(object$coefficients))<=treat_pos)>0){
    num_NA = sum(which(is.na(object$coefficients))<=treat_pos) # number of NA coefficients before TE -> del. cols before TE
    C[treat_pos-num_NA] = 1 # if some cols were deleted before the TE, we need a different row of the resp. matrix in the end
  } else {stop("Check for multicollinearity")} # just in case something is wrong -> but shouldn't be the case -> can be deleted after some testing
  
  # Access qr decomposition
  Q = qr.Q(object$qr)
  R = qr.R(object$qr)
  n_NA = sum(is.na(object$coefficients)) # we need to leave the last columns and rows out to get the OLS coefficients if there are NAs
  ncol = dim(R)[1]
  
  #### Calculate the outcome weights:
  omega = C %*% backsolve(R[1:(ncol-n_NA),1:(ncol-n_NA)],t(Q[,1:(ncol-n_NA)]))# this is a numerically stable way and follows how the lm command computes the inverse of X'X. For an intuitive example
  # on how numerical instabilities can arise check e.g., https://genomicsclass.github.io/book/pages/qr_and_regression.html
  # The last rows/cols are left out if there is multicollinearity, i.e. the matrix X'X is singular
  # This is how we would get back the OLS coefficients (by multiplying with Y)
  
  # if the original command uses sample weights we need to multiply by sqrt(w) (elementwise) (the omegas)
  # (as weighted lm multiplies X and Y by sqrt(w))
  # additionally, if some observations have weight 0 they are "deleted" for the lm estimation.
  # We need to take care of this to make the dimensions match
  if(is.numeric(object$model$`(weights)`)){
    w = object$model$`(weights)`
    ok = w!=0
    w = w[ok]
    omega_help = omega*sqrt(w) # after doing this we can multiply by Y to get back the coefficient
    omega = rep(0,n)
    omega[ok] = omega_help
  }
  
  Y = as.matrix(model.frame(object)[,1],n)
  if(!isTRUE(all.equal(as.numeric(omega %*% Y),
                       as.numeric(object$coefficients[treat_pos])))){
    warning("Estimated Treatment Effects using weights differ from original estimates.")
  }
  
  output = list(
    "omega" = matrix(omega,nrow=1),
    "treat" = D
  )
  
  class(output) = c("get_outcome_weights")
  
  return(output)
}



#' Outcome weights for ivreg
#'
#' @description Post-estimation command to calculate the outcome weights for
#' Instrumental-Variable Regression implemented via the \code{\link[AER]{ivreg}} command.
#' Automatically takes sampling weights into account if used.
#'
#' @param object An object of class \code{ivreg}, i.e. the result of running \code{\link[AER]{ivreg}}.
#'
#' @param treat_pos The position of the treatment indicator in the model matrix. By default 2,
#' which is appropriate for the case with an intercept and the (endogneous) treatment indicator as the first explanatory variable.
#' @param ... Pass potentially generic \link{get_outcome_weights} options.
#'
#' @return A vector of implied weights.
#'
#' @export
#'
get_outcome_weights.ivreg = function(object,treat_pos = 2,...){
  Y = as.matrix(model.frame(object)[,1])

  X_const = model.matrix(object$terms$regressors, object$model, contrasts = object$contrasts$regressors)
  Z_const = model.matrix(object$terms$instruments, object$model, contrasts = object$contrasts$instruments)
  
  C = matrix(0,1,ncol(X_const))
  C[treat_pos] = 1
  D = X_const[,treat_pos]
  
  # sampling weight if necessary
  if (!is.null(object$model$`(weights)`)&is.numeric(object$model$`(weights)`)){
    w = object$model$`(weights)`
    X_const = X_const*sqrt(w)
    Z_const = Z_const*sqrt(w)
    omega = C%*%solve(t(X_const)%*%Z_const%*%solve(t(Z_const)%*%Z_const)%*%t(Z_const)%*%X_const)%*%t(X_const)%*%Z_const%*%solve(t(Z_const)%*%Z_const)%*%t(Z_const)*sqrt(w)
    
  } else{
    omega = C%*%solve(t(X_const)%*%Z_const%*%solve(t(Z_const)%*%Z_const)%*%t(Z_const)%*%X_const)%*%t(X_const)%*%Z_const%*%solve(t(Z_const)%*%Z_const)%*%t(Z_const)
  }
  
  if(!isTRUE(all.equal(as.numeric(omega %*% Y),
                       as.numeric(object$coefficients[treat_pos])))){
    warning("Estimated Treatment Effects using weights differ from original estimates.")
  }
  
  output = list(
    "omega" = matrix(omega,nrow=1),
    "treat" = D
  )
  
  class(output) = c("get_outcome_weights")
  
  return(output)
}
