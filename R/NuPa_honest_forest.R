#' Nuisance parameter estimation via honest random forest
#' 
#' This function estimates different nuisance parameters using the honest random forest
#' implementation of the 'grf' package
#'
#' @param NuPa String vector specifying the nuisance parameters to be estimated.
#' Currently supported: \code{c("Y.hat","Y.hat.d","Y.hat.z","D.hat","D.hat.z","Z.hat")} 
#' @param X Covariate matrix with N rows and p columns.
#' @param Y Optional numeric vector containing the outcome variable.
#' @param D Optional binary treatment variable.
#' @param Z Optional binary instrumental variable.
#' @param n_cf_folds Number of cross-fitting folds. Default is 5.
#' @param n_reps Number of repetitions of cross-fitting. Default is 1.
#' @param cluster Optional vector of cluster variable if cross-fitting should account for clusters.
#' @param progress If TRUE, progress of nuisance parameter estimation reported.
#' @param ... Options passed to the \code{\link[grf]{regression_forest}}.
#'
#' @return List of two lists. 
#' - \code{predictions} contains the requested nuisance parameters
#' - \code{smoothers} contains the smoother matrices of requested outcome nuisance parameters
#' - \code{cf_mat} Array of dimension n_reps x N x n_cf_folds storing indicators representing the folds used in estimation.
#' 
#' @references 
#' Wager, S., & Athey, S. (2018). Estimation and inference of heterogeneous treatment effects using random forests. 
#' Journal of the American Statistical Association, 113(523), 1228-1242.
#'
#' @export
#'
NuPa_honest_forest = function(NuPa = c("Y.hat","Y.hat.d","Y.hat.z","D.hat","D.hat.z","Z.hat"),
                              X, 
                              Y=NULL, 
                              D=NULL, 
                              Z=NULL,
                              n_cf_folds=5,
                              n_reps=1,
                              cluster=NULL,
                              progress=FALSE,
                              ...) {
  # Sanity checks
  supported_NuPas = c("Y.hat","Y.hat.d","Y.hat.z","D.hat","D.hat.z","Z.hat")
  not_supported = NuPa[!NuPa %in% supported_NuPas]
  if (length(not_supported) > 0) {
    stop(paste("Error: The following nuisance parameters specified in NuPa are not supported:", 
               paste(not_supported, collapse = ", ")))}
  if (is.null(Y) & any(NuPa %in% c("Y.hat","Y.hat.d","Y.hat.z"))) {
    stop("Please specify Y if at least one of c(\"Y.hat\",\"Y.hat.d\",\"Y.hat.z\") is specified in NuPa")}
  if (is.null(D) & any(NuPa %in% c("Y.hat.d","D.hat","D.hat.z"))) {
    stop("Please specify D if at least one of c(\"Y.hat.d\",\"D.hat\",\"D.hat.z\") is specified in NuPa")}
  if (is.null(Z) & any(NuPa %in% c("Y.hat.z","D.hat.z","Z.hat"))) {
    stop("Please specify Z if at least one of c(\"Y.hat.z\",\"D.hat.z\",\"Z.hat\") is specified in NuPa")}
  
  # Preps
  N = nrow(X)
  
  # Initialize nuisance parameters for all
  Y.hat = Y.hat.d0 = Y.hat.d1 = Y.hat.z0 = Y.hat.z1 = 
    D.hat = D.hat.z0 = D.hat.z1 = Z.hat = 
    S = S.d0 = S.d1 = S.z0 = S.z1 = 
    "This nuisance parameter was not specified and is therefore empty."
  
  # Initialize nuisance parameters and smoother matrices to be filled
  if ("Y.hat" %in% NuPa) {Y.hat = matrix(NA,N,n_reps); S = array(0,c(n_reps,N,N))}
  if ("Y.hat.d" %in% NuPa) {Y.hat.d0 = Y.hat.d1 = matrix(NA,N,n_reps); S.d0 = S.d1 = array(0,c(n_reps,N,N))}
  if ("Y.hat.z" %in% NuPa) {Y.hat.z0 = Y.hat.z1 = matrix(NA,N,n_reps); S.z0 = S.z1 = array(0,c(n_reps,N,N))}
  if ("D.hat" %in% NuPa) D.hat = matrix(NA,N,n_reps)
  if ("D.hat.z" %in% NuPa) D.hat.z0 = D.hat.z1 = matrix(NA,N,n_reps)
  if ("Z.hat" %in% NuPa) Z.hat = matrix(NA,N,n_reps)
  cf_mat = array(NA,c(n_reps,N,n_cf_folds))
  
  for (r in 1:n_reps) {
    if (isTRUE(progress)) cat(sprintf("[%s]: repetition %s started\n", 
                                      format(Sys.time(), "%Y-%m-%d %H:%M:%S"), toString(r)))
    # Get cross-fitting folds
    cf_mat[r,,] = prep_cf_mat(N,n_cf_folds,cl=cluster) 
    
    # Gradually fill empty vectors and matrices
    for (i in 1:n_cf_folds){
      if (isTRUE(progress)) cat(sprintf("%s [%s]: cross-fitting fold %s of repetition %s started\n", 
                  "    ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), toString(i), toString(r)))
      train_fold = !cf_mat[r,,i]
      pred_fold = cf_mat[r,,i]
      if ("Y.hat" %in% NuPa) {
        if (isTRUE(progress)) cat(sprintf("%s [%s]: estimation of NuPa %s started\n", 
                                          "         ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "Y.hat"))
        rf_Yhat = regression_forest(X[train_fold,], Y[train_fold], compute.oob.predictions = F, ...)
        Y.hat[pred_fold,r] = predict(rf_Yhat, X[pred_fold,])$predictions
        S[r,pred_fold, train_fold] = as.matrix(grf::get_forest_weights(rf_Yhat, X[pred_fold,]))
      }
      if ("Y.hat.d" %in% NuPa) {
        if (isTRUE(progress)) cat(sprintf("%s [%s]: estimation of NuPa %s started\n", 
                                          "         ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "Y.hat.d"))
        rf_Yhatd0 = regression_forest(X[train_fold & D==0,], Y[train_fold & D==0], compute.oob.predictions = F, ...)
        Y.hat.d0[pred_fold,r] = predict(rf_Yhatd0, X[pred_fold,])$predictions
        S.d0[r,pred_fold, train_fold & D==0] = as.matrix(grf::get_forest_weights(rf_Yhatd0, X[pred_fold,]))
        
        rf_Yhatd1 = regression_forest(X[train_fold & D==1,], Y[train_fold & D==1], compute.oob.predictions = F, ...)
        Y.hat.d1[pred_fold,r] = predict(rf_Yhatd1, X[pred_fold,])$predictions
        S.d1[r,pred_fold, train_fold  & D==1] = as.matrix(grf::get_forest_weights(rf_Yhatd1, X[pred_fold,]))
      }
      if ("Y.hat.z" %in% NuPa) {
        if (isTRUE(progress)) cat(sprintf("%s [%s]: estimation of NuPa %s started\n", 
                                          "         ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "Y.hat.z"))
        rf_Yhatz0 = regression_forest(X[train_fold & Z==0,], Y[train_fold & Z==0], compute.oob.predictions = F, ...)
        Y.hat.z0[pred_fold,r] = predict(rf_Yhatz0, X[pred_fold,])$predictions
        S.z0[r,pred_fold, train_fold  & Z==0] = as.matrix(grf::get_forest_weights(rf_Yhatz0, X[pred_fold,]))
        
        rf_Yhatz1 = regression_forest(X[train_fold & Z==1,], Y[train_fold & Z==1], compute.oob.predictions = F, ...)
        Y.hat.z1[pred_fold,r] = predict(rf_Yhatz1, X[pred_fold,])$predictions
        S.z1[r,pred_fold, train_fold  & Z==1] = as.matrix(grf::get_forest_weights(rf_Yhatz1, X[pred_fold,]))
      }
      if ("D.hat" %in% NuPa) {
        if (isTRUE(progress)) cat(sprintf("%s [%s]: estimation of NuPa %s started\n", 
                                          "         ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "D.hat"))
        rf_Dhat = regression_forest(X[train_fold,], D[train_fold], compute.oob.predictions = F, ...)
        D.hat[pred_fold,r] = predict(rf_Dhat, X[pred_fold,])$predictions
      }
      if ("D.hat.z" %in% NuPa) {
        if (isTRUE(progress)) cat(sprintf("%s [%s]: estimation of NuPa %s started\n", 
                                          "         ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "D.hat.z"))
        rf_Dhatz0 = regression_forest(X[train_fold & Z==0,], D[train_fold & Z==0], compute.oob.predictions = F, ...)
        D.hat.z0[pred_fold,r] = predict(rf_Dhatz0, X[pred_fold,])$predictions
        
        rf_Dhatz1 = regression_forest(X[train_fold & Z==1,], D[train_fold & Z==1], compute.oob.predictions = F, ...)
        D.hat.z1[pred_fold,r] = predict(rf_Dhatz1, X[pred_fold,])$predictions
      }
      if ("Z.hat" %in% NuPa) {
        if (isTRUE(progress)) cat(sprintf("%s [%s]: estimation of NuPa %s started\n", 
                                          "         ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "Z.hat"))
        rf_Zhat = regression_forest(X[train_fold,], Z[train_fold], compute.oob.predictions = F, ...)
        Z.hat[pred_fold,r] = predict(rf_Zhat, X[pred_fold,])$predictions
      }
    } # end loop over i
    
  }
  
  list_predictions = list("Y.hat"=Y.hat,"Y.hat.d0"=Y.hat.d0,"Y.hat.d1"=Y.hat.d1,
                          "Y.hat.z0"=Y.hat.z0,"Y.hat.z1"=Y.hat.z1,"D.hat"=D.hat,
                          "D.hat.z0"=D.hat.z0,"D.hat.z1"=D.hat.z1,"Z.hat"=Z.hat)
  list_smoothers = list("S"=S,"S.d0"=S.d0, "S.d1"=S.d1, "S.z0"=S.z0, "S.z1"=S.z1)
  
  list("predictions"=list_predictions, "smoothers"=list_smoothers, "cf_mat" = cf_mat)
}

