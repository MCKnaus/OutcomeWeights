######### causal forest #########

#' Outcome weights for the \code{\link[grf]{causal_forest}} function
#'
#' @description Post-estimation command to extract outcome weights for causal forest
#' implemented via the \code{\link[grf]{causal_forest}} function from the \pkg{grf} package.
#'
#' @param object An object of class \code{causal_forest}, i.e. the result of running \code{\link[grf]{causal_forest}}.
#' @param ... Pass potentially generic \link{get_outcome_weights} options.
#' @param S A smoother matrix reproducing the outcome predictions used in building the \code{\link[grf]{instrumental_forest}}. 
#' Obtained by calling \code{get_forest_weights()} for the \code{\link[grf]{regression_forest}} object producing the outcome predictions.
#' @param newdata Corresponds to \code{newdata} option in \code{\link[grf]{predict.causal_forest}}. If \code{NULL}, 
#' out-of-bag outcome weights, otherwise for those for the provided test data returned.
#' @param S.tau Required if \code{target != "CATE"}, then S.tau is the CATE smoother obtained from running \code{get_outcome_weights()}
#' with \code{target == "CATE"}.
#' @param target Target parameter for which outcome weights should be extracted. Currently \code{c("CATE","ATE")} implemented.
#' @param checks Default \code{TRUE} checks whether weights numerically replicate original estimates. Only set \code{FALSE} if you 
#' know what you are doing and need to save computation time.
#'
#' @return \link{get_outcome_weights} object with `omega` containing weights and `treat` the treatment
#' 
#' @examples
#' \donttest{
#' # Sample from DGP borrowed from grf documentation
#' n = 500
#' p = 10
#' X = matrix(rnorm(n * p), n, p)
#' W = rbinom(n, 1, 0.5)
#' Y = pmax(X[, 1], 0) * W + X[, 2] + pmin(X[, 3], 0) + rnorm(n)
#' 
#' # Run outcome regression and extract smoother matrix
#' forest.Y = grf::regression_forest(X, Y)
#' Y.hat = predict(forest.Y)$predictions
#' outcome_smoother = grf::get_forest_weights(forest.Y)
#' 
#' # Run causal forest with external Y.hats
#' c.forest = grf::causal_forest(X, Y, W, Y.hat = Y.hat)
#'
#' # Predict on out-of-bag training samples.
#' cate.oob = predict(c.forest)$predictions
#' 
#' # Predict using the forest.
#' X.test = matrix(0, 101, p)
#' X.test[, 1] = seq(-2, 2, length.out = 101)
#' cate.test = predict(c.forest, X.test)$predictions
#' 
#' # Calculate outcome weights
#' omega_oob = get_outcome_weights(c.forest,S = outcome_smoother)
#' omega_test = get_outcome_weights(c.forest,S = outcome_smoother,newdata = X.test)
#' 
#' # Observe that they perfectly replicate the original CATEs
#' all.equal(as.numeric(omega_oob$omega %*% Y), 
#'           as.numeric(cate.oob))
#' all.equal(as.numeric(omega_test$omega %*% Y), 
#'           as.numeric(cate.test))
#' 
#' # Also the ATE estimates are perfectly replicated
#' omega_ate = get_outcome_weights(c.forest,target = "ATE", 
#'                                 S = outcome_smoother, 
#'                                 S.tau = omega_oob$omega)
#' all.equal(as.numeric(omega_ate$omega %*% Y),
#'           as.numeric(grf::average_treatment_effect(c.forest, target.sample = "all")[1]))
#' 
#' # The omega weights can be plugged into balancing packages like cobalt
#' }
#' 
#' @references 
#' Athey, S., Tibshirani, J., & Wager, S. (2019). Generalized random forest. The Annals of Statistics, 47(2), 1148-1178.
#' 
#' Knaus, M. C. (2024). Treatment effect estimators as weighted outcomes, \url{https://arxiv.org/abs/2411.11559}.
#' 
#' @importFrom stats predict
#' @importFrom methods as
#' @importFrom grf regression_forest
#'
#' @export
#' @method get_outcome_weights causal_forest
#' 
get_outcome_weights.causal_forest = function(object,...,
                                             S,
                                             newdata=NULL,
                                             S.tau=NULL,
                                             target="CATE",
                                             checks=TRUE){
  ### Extract and define important components
  D = object$W.orig
  Y = object$Y.orig
  Dhat = object$W.hat
  Dres = D - Dhat
  n = length(Y)
  ones = matrix(1,n,1)
  
    ### Sanity checks
  if (!(target %in% c("ATE", "CATE"))) stop("Currently online CATE and ATE as target parameters available.")
  if (is.null(S.tau) & target == "ATE") stop("You specify target == \"ATE\" but do not provide S.tau as the CATE smoother matrix.
                                            Please run get_outcome_weights first with target == \"CATE\" and pass the results as S.tau.")
  if (!is.null(newdata) & target == "ATE") warning("newdata ignored when calculating ATE weights.")
  if (checks) {
    Y.hat_weights = S %*% Y
    if (!all.equal(as.numeric(object$Y.hat), as.numeric(Y.hat_weights))) stop("Outcome smoother matrix does not replicate Y.hat used by causal forest.")
  }
  
  if (target == "CATE") {
    # make the matrix a sparse matrix such that the C++ function accepts it
    if (!inherits(S, "dgCMatrix")) S = as(S,"dgCMatrix") 
    
    ### Get the causal forest CATE weights
    if (is.null(newdata)) alpha = grf::get_forest_weights(object)
    else alpha = grf::get_forest_weights(object, newdata)
    
    # call the C++ function that calculates the Ztildex matrix:
    scaled_Ztildex = scaled_Ztildex_maker(alpha,Dres,Dres)  # special case of instrumental forest with Zres = Dres0
    # Transformation matrix
    Tcf = (diag(n) - S)
    # Calculate the omega matrix
    omega = as.matrix(scaled_Ztildex %*% Tcf)
  
    if (checks) {
      if (is.null(newdata)) cates_cf = predict(object)$predictions
      else cates_cf = predict(object,newdata=newdata)$predictions
      if(!isTRUE(all.equal(as.numeric(omega %*% Y),as.numeric(cates_cf)))){
        warning("CATEs produced using weights differ from original estimates.")
      }
    }
  }  # end if (target == "CATE")
  
  
  if (target == "ATE") {
    lambda1 = D / Dhat
    lambda0 = (1-D) / (1-Dhat)

    # Element-wise operations for scaling
    scaled_S.tau = (D - Dhat) * S.tau  # Element-wise multiplication

    # Calculate S adjustment
    S_adjusted = diag(n) - S - scaled_S.tau

    # Compute the final matrix
    omage_aipw_grf = S.tau + (lambda1 - lambda0) * S_adjusted  # Element-wise multiplication

    omega = matrix(t(ones) %*% omage_aipw_grf/ n,nrow=1)
    
    if (checks) {
      cates_cf = predict(object)$predictions
      if(!isTRUE(all.equal(as.numeric(omega %*% Y),
                           as.numeric(grf::average_treatment_effect(object, target.sample = "all")[1])))){
        warning("Estimated Treatment Effects using weights differ from original estimates.")
      }
    }
  } # end if (target == "ATE")
  
  output = list(
    "omega" = omega,
    "treat" = D
  )
  
  class(output) = c("get_outcome_weights")
  
  return(output)
}


######### instrumental forest #########

#' Outcome weights for the \code{\link[grf]{instrumental_forest}} function
#'
#' @description Post-estimation command to extract outcome weights for instrumental forest
#' implemented via the \code{\link[grf]{instrumental_forest}} function from the \pkg{grf} package.
#' 
#' @param object An object of class \code{instrumental_forest}, i.e. the result of running \code{\link[grf]{instrumental_forest}}.
#' @param ... Pass potentially generic \link{get_outcome_weights} options.
#' @param S A smoother matrix reproducing the outcome predictions used in building the \code{\link[grf]{instrumental_forest}}. 
#' Obtained by calling \code{get_forest_weights()} for the \code{\link[grf]{regression_forest}} object producing the outcome predictions.
#' @param newdata Corresponds to \code{newdata} option in \code{\link[grf]{predict.instrumental_forest}}. If \code{NULL}, 
#' out-of-bag outcome weights, otherwise for those for the provided test data returned.
#' @param checks Default \code{TRUE} checks whether weights numerically replicate original estimates. Only set \code{FALSE} if you 
#' know what you are doing and want to save computation time.
#'
#' @return \link{get_outcome_weights} object with `omega` containing weights and `treat` the treatment
#' 
#' @examples
#' \donttest{
#' # Sample from DGP borrowed from grf documentation
#' n = 2000
#' p = 5
#' X = matrix(rbinom(n * p, 1, 0.5), n, p)
#' Z = rbinom(n, 1, 0.5)
#' Q = rbinom(n, 1, 0.5)
#' W = Q * Z
#' tau =  X[, 1] / 2
#' Y = rowSums(X[, 1:3]) + tau * W + Q + rnorm(n)
#' 
#' # Run outcome regression and extract smoother matrix
#' forest.Y = grf::regression_forest(X, Y)
#' Y.hat = predict(forest.Y)$predictions
#' outcome_smoother = grf::get_forest_weights(forest.Y)
#' 
#' # Run instrumental forest with external Y.hats
#' iv.forest = grf::instrumental_forest(X, Y, W, Z, Y.hat = Y.hat)
#'
#' # Predict on out-of-bag training samples.
#' iv.pred = predict(iv.forest)$predictions
#' 
#' omega_if = get_outcome_weights(iv.forest, S = outcome_smoother)
#' 
#' # Observe that they perfectly replicate the original CLATEs
#' all.equal(as.numeric(omega_if$omega %*% Y), 
#'           as.numeric(iv.pred))
#'
#' }
#' 
#' @references 
#' Athey, S., Tibshirani, J., & Wager, S. (2019). Generalized random forest. The Annals of Statistics, 47(2), 1148-1178.
#' 
#' Knaus, M. C. (2024). Treatment effect estimators as weighted outcomes, \url{https://arxiv.org/abs/2411.11559}.
#' 
#' @importFrom stats predict
#' @importFrom methods as
#' @importFrom grf regression_forest
#'
#' @export
#' @method get_outcome_weights instrumental_forest
#' 
get_outcome_weights.instrumental_forest = function(object,...,
                                                   S,
                                                   newdata=NULL,
                                                   checks=TRUE){
  ### Extract and define important components
  D = object$W.orig
  Y = object$Y.orig
  Dres = D - object$W.hat
  Zres = object$Z.orig - object$Z.hat
  n = length(Y)
  ones = matrix(1,n,1)
  
  ### Sanity checks
  if (checks) {
    Y.hat_weights = S %*% Y
    if (!all.equal(as.numeric(object$Y.hat), as.numeric(Y.hat_weights))) stop("Outcome smoother matrix does not replicate Y.hat used by causal forest.")
  }
  
  # make the matrix a sparse matrix such that the C++ function accepts it
  if (!inherits(S, "dgCMatrix")) S = as(S,"dgCMatrix") 
  
  ### Get the causal forest CATE weights
  if (is.null(newdata)) alpha = grf::get_forest_weights(object)
  else alpha = grf::get_forest_weights(object, newdata)
  
  # call the C++ function that calculates the Ztildex matrix:
  scaled_Ztildex = scaled_Ztildex_maker(alpha,Zres,Dres)  # special case of instrumental forest with Zres = Dres0
  # Transformation matrix
  Tcf = (diag(n) - S)
  # Calculate the omega matrix
  omega = as.matrix(scaled_Ztildex %*% Tcf)
  
  if (checks) {
    if (is.null(newdata)) clates_cf = predict(object)$predictions
    else clates_cf = predict(object,newdata=newdata)$predictions
    if(!isTRUE(all.equal(as.numeric(omega %*% Y),as.numeric(clates_cf)))){
      warning("CATEs produced using weights differ from original estimates.")
    }
  }
  
  output = list(
    "omega" = omega,
    "treat" = D
  )
  
  class(output) = c("get_outcome_weights")
  
  return(output)
}


