#' Outcome weights for a \code{DoubleML} object
#'
#' @description Post-estimation command to extract outcome weights for a 
#' \code{DoubleML} object from the \pkg{DoubleML} package.
#'
#' @param object A \code{DoubleML} object of class \code{DoubleMLPLR}, 
#'               \code{DoubleMLPLIV}, \code{DoubleMLIIVM}, or \code{DoubleMLIRM}.
#' @param dml_data The \code{DoubleMLData} object used in estimation.
#' @param ... Pass potentially generic \link{get_outcome_weights} options.
#'
#' @return 
#' An object of class \code{get_outcome_weights} containing:
#' \describe{
#'   \item{omega}{Numeric vector of outcome weights of length \eqn{N}.}
#'   \item{treat}{Numeric vector of treatment assignments.}
#' }
#' 
#' @details
#' Outcome weights are available for \code{DoubleMLPLR}, \code{DoubleMLPLIV},
#' \code{DoubleMLIIVM}, and \code{DoubleMLIRM} models.
#' 
#' @examples
#' \dontrun{
#' set.seed(123)
#'
#' # Set the parameters
#' params <- list(alpha = 0, subsample = 1, max_delta_step = 0, base_score = 0)
#'
#' # Define DML objects
#' ml_regr_xgb <- do.call(lrn, c(list("regr.xgboost"), params))
#' ml_l_plr <- ml_regr_xgb$clone()
#' ml_m_plr <- ml_regr_xgb$clone()
#' plr_data <- DoubleML::make_plr_CCDDHNR2018(500, dim_x = 10)
#' plr_obj <- DoubleML::DoubleMLPLR$new(plr_data, ml_l_plr, ml_m_plr, n_folds = 3)
#'
#' # Fit the model
#' plr_obj$fit(store_models = TRUE, store_predictions = TRUE)
#'
#' # Calculate the outcome weights
#' omega_plr <- get_outcome_weights(plr_obj, plr_data)
#' 
#' # Check equivalence 
#' all.equal(omega_plr$omega %*% plr_data$data$y, plr_obj$all_coef)
#' }
#'
#' @references 
#' Knaus, M. C. (2024). Treatment effect estimators as weighted outcomes, \url{https://arxiv.org/abs/2411.11559}.
#' 
#' Bach P, Kurz MS, Chernozhukov V, Spindler M, Klaassen S (2024). “DoubleML: 
#' An Object-Oriented Implementation of Double Machine Learning in R.” 
#' Journal of Statistical Software, 108(3), 1-56, \url{https://doi.org/10.18637/jss.v108.i03}.
#'
#' @export
get_outcome_weights.DoubleML = function(object, dml_data, ...) {
  ## Preps
  type  <- class(object)[1]
  data  <- as.matrix(dml_data$data)
  preds <- lapply(object$predictions, as.numeric)
  
  y_col  <- dml_data$y_col
  x_cols <- dml_data$x_cols
  d_col  <- dml_data$d_cols
  z_col  <- dml_data$z_cols
  
  Y <- data[, y_col]
  X <- data[, x_cols, drop = FALSE]
  D <- if (!is.null(d_col)) data[, d_col]
  Z <- if (!is.null(z_col)) data[, z_col]
  N <- nrow(data)
  
  test_ids_list  <- object$smpls[[1]]$test_ids
  train_ids_list <- object$smpls[[1]]$train_ids
  
  
  ## Get outcome weights
  ### PLR ###
  if (type == "DoubleMLPLR") {
    fold_models <- object$models$ml_l$d[[1]]
    
    S <- get_DoubleML_smoother(
      fold_models, test_ids_list, train_ids_list, data, x_cols, y_col, NULL
    )
    
    Z.tilde = D.tilde = D - preds$ml_m
    T_mat = diag(N) - S
    omega = pive_weight_maker(Z.tilde, D.tilde, T_mat)
    
    
    ### PLR-IV ###
  } else if (type == "DoubleMLPLIV") {
    fold_models <- object$models$ml_l$d[[1]]
    
    S <- get_DoubleML_smoother(
      fold_models, test_ids_list, train_ids_list, data, x_cols, y_col, NULL
    )
    
    Z.tilde = Z - preds$ml_m
    D.tilde = D - preds$ml_r
    T_mat = diag(N) - S
    omega = pive_weight_maker(Z.tilde, D.tilde, T_mat)
    
    
    ### AIPW-ATE ###
  } else if (type == "DoubleMLIRM") {
    fold_models_g0 <- object$models$ml_g0$d[[1]]
    fold_models_g1 <- object$models$ml_g1$d[[1]]
    
    sub_d0 <- D == 0
    sub_d1 <- D == 1
    
    S.d0 <- get_DoubleML_smoother(
      fold_models_g0, test_ids_list, train_ids_list, data, x_cols, y_col, sub_d0
    )
    S.d1 <- get_DoubleML_smoother(
      fold_models_g1, test_ids_list, train_ids_list, data, x_cols, y_col, sub_d1
    )
    
    Z.tilde = D.tilde = rep(1, N)
    lambda1 = D / preds$ml_m
    lambda0 = (1 - D) / (1 - preds$ml_m)
    T_mat = S.d1 - S.d0 + 
      lambda1 * (diag(N) - S.d1) - 
      lambda0 * (diag(N) - S.d0)
    omega = pive_weight_maker(Z.tilde, D.tilde, T_mat)
    
    
    ### Wald-AIPW ###
  } else if (type == "DoubleMLIIVM") {
    fold_models_g0 <- object$models$ml_g0$d[[1]]
    fold_models_g1 <- object$models$ml_g1$d[[1]]
    
    sub_z0 <- Z == 0
    sub_z1 <- Z == 1
    
    S.z0 <- get_DoubleML_smoother(
      fold_models_g0, test_ids_list, train_ids_list, data, x_cols, y_col, sub_z0
    )
    S.z1 <- get_DoubleML_smoother(
      fold_models_g1, test_ids_list, train_ids_list, data, x_cols, y_col, sub_z1
    )
    
    lambdaz1 = Z / preds$ml_m
    lambdaz0 = (1 - Z) / (1 - preds$ml_m)
    Z.tilde = rep(1, N)
    D.tilde = preds$ml_r1 - preds$ml_r0 + 
      lambdaz1 * (D - preds$ml_r1) - 
      lambdaz0 * (D - preds$ml_r0)
    T_mat = S.z1 - S.z0 + 
      lambdaz1 * (diag(N) - S.z1) - 
      lambdaz0 * (diag(N) - S.z0)
    omega = pive_weight_maker(Z.tilde, D.tilde, T_mat)
    
  } else {
    stop(
      "Outcome weights are extracted for DoubleMLPLR, DoubleMLPLIV, DoubleMLIIVM, and DoubleMLIRM models"
    )
  }
  
  output = list("omega" = omega, "treat" = D)
  class(output) = c("get_outcome_weights")
  return(output)
}


#' Create smoother matrix for \code{DoubleML} models
#'
#' Computes an \eqn{N \times N} smoother matrix for the corresponding nuisance parameter
#' by looping over the fold-specific models trained within \code{DoubleML}.
#'
#' @param fold_models List of fold-specific models from the \code{DoubleML} object.
#' @param test_ids_list List of test indices for each fold.
#' @param train_ids_list List of training indices for each fold.
#' @param data A numeric matrix with the variables used to construct
#'             the \code{DoubleMLData} object.
#' @param subset Logical vector indicating which observations to use for
#'  extracting smoother weights. If not provided, all observations are used.
#' @param x_cols Character vector of covariate column names.
#' @param y_col Character string specifying the outcome column name.
#'
#' @return An \eqn{N \times N} smoother matrix.
#' 
#' @references 
#' Knaus, M. C. (2024). Treatment effect estimators as weighted outcomes, \url{https://arxiv.org/abs/2411.11559}.
#' 
#' Bach P, Kurz MS, Chernozhukov V, Spindler M, Klaassen S (2024). “DoubleML: 
#' An Object-Oriented Implementation of Double Machine Learning in R.” 
#' Journal of Statistical Software, 108(3), 1-56, \url{https://doi.org/10.18637/jss.v108.i03}.
#'
#' @keywords internal
#'
get_DoubleML_smoother <- function(fold_models,
                                  test_ids_list,
                                  train_ids_list,
                                  data,
                                  x_cols,
                                  y_col,
                                  subset = NULL) {
  N <- nrow(data)
  S <- matrix(0, N, N)
  
  if (is.null(subset)) subset <- rep(TRUE, N)
  
  for (i in seq_along(fold_models)) {
    model <- fold_models[[i]]$model
    
    # ranger's double-wrapping quirk: 
    if (inherits(model$model, "ranger")) model <- model$model
    
    # xgboost's feature ordering quirk: 
    if (inherits(model, "xgb.Booster")) x_cols <- xgboost::getinfo(model, "feature_name")
    
    test_ids  <- test_ids_list[[i]]
    train_ids <- train_ids_list[[i]]
    train_ids_sub <- train_ids[subset[train_ids]]
    
    test_data   <- data[test_ids,      x_cols, drop = FALSE]
    train_data  <- data[train_ids_sub, x_cols, drop = FALSE]
    train_label <- data[train_ids_sub, y_col,  drop = FALSE]
    
    S_fold <- get_smoother_weights(model, X = train_data, Y = train_label, Xnew = test_data)
    
    S[test_ids, train_ids_sub] <- S_fold
  }
  
  return(S)
}
