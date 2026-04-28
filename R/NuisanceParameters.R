#' Smoother weights for outcome models
#'
#' Extract smoother matrices for nuisance parameter outcome models from a
#' \code{NuisanceParameters} object created by \code{\link[NuisanceParameters]{nuisance_parameters}}.
#'
#' @param object A \code{NuisanceParameters} object created by
#'   \code{\link[NuisanceParameters]{nuisance_parameters}}.
#' @param NuPa Character vector specifying which nuisance parameters to extract
#'  smoothers for. Available outcome nuisance parameters
#'  are \code{c("Y.hat", "Y.hat.d", "Y.hat.z")}.
#' @param quiet If \code{FALSE}, the method currently computed is printed to the console.
#' @param subset Logical vector indicating which observations to use when determining
#'   ensemble weights. If not provided, all observations are used.
#' @param ... Additional arguments.
#'
#' @return A list of matrices of dimension \code{nrow(X)} × \code{nrow(X)} containing ensemble smoother
#'   matrices
#'
#' @method get_smoother_weights NuisanceParameters
#' @export
get_smoother_weights.NuisanceParameters <- function(object,
                                                    NuPa,
                                                    quiet = TRUE,
                                                    subset = NULL,
                                                    ...) {
  # Input checks
  NuPa <- match.arg(NuPa, choices = c("Y.hat", "Y.hat.d", "Y.hat.z"), several.ok = TRUE)
  
  if (!inherits(object, "NuisanceParameters")) {
    stop("object must be a NuisanceParameters object (i.e., the output of nuisance_parameters())")
  }
  
  X <- object$numbers$X
  Y <- object$numbers$Y
  cv <- object$numbers$cv
  d_mat <- object$numbers$d_mat
  z_mat <- object$numbers$z_mat
  cf_mat <- object$numbers$cf_mat
  models <- object$models
  
  methods <- object$numbers$methods
  methods <- methods[names(methods) %in% NuPa]
  
  all_methods  <- unique(unlist(lapply(methods, function(m) unlist(lapply(m, `[[`, "method")))))
  smooth_support <- c("mean", "ols", "ridge", "plasso", "forest_grf", 
                      "knn", "forest_drf", "xgboost", "rlasso", "ranger")
  
  no_smooth <- setdiff(all_methods, smooth_support)
  if (length(no_smooth) > 0) {
    stop(paste("Smoother matrices cannot be extracted for method(s):", 
               paste(no_smooth, collapse = ", "),
               "- refer to documentation for supported methods."))
  }
  
  # XGBoost check
  if ("xgboost" %in% all_methods) {
    message(
      "Note that numerical equivalence is only at single precision, not at double precision",
      " because XGBoost contributes to the predictions. This is usually no reason to worry."
    )
  }
  
  # Models check
  m_names <- paste0(NuPa, "_m")
  outcome_m <- m_names[m_names %in% names(models) & !vapply(models[m_names], is.character, logical(1))]
  missing_m <- setdiff(m_names, outcome_m)
  
  if (length(missing_m) > 0) {
    message(
      "Missing models: ", paste(missing_m, collapse = ", "),
      "\nContinuing with outcome models: ", paste(outcome_m, collapse = ", ")
    )
  }
  if (length(outcome_m) == 0) stop("No valid (ensemble) models found in object$models.")
  if (is.null(subset)) subset <- rep(TRUE, length(Y))
  
  results <- list()
  
  # Process each model
  for (m_name in outcome_m) {
    current_model <- models[[m_name]]
    nupa <- sub("_m$", "", m_name)
    nupa_smoother <- list()
    
    ## Short stacking
    # Case 1: Single ensemble model (i.e. Y.hat)
    if (cv == 1) {
      if ("ens_object" %in% names(current_model)) {
        nupa_smoother <- get_smoother(
          object = current_model,
          methods = methods[[nupa]],
          X = X,
          cv = cv,
          cf_mat = cf_mat,
          quiet = quiet
        )
        # Case 2: List of models (i.e. Y.hat.d, Y.hat.z)
      } else {
        for (i in seq_along(current_model)) {
          sub_element <- current_model[[i]]
          subset <- get_subset(sub_element, models, d_mat, z_mat)
          
          nupa_smoother[[i]] <- get_smoother(
            object = sub_element,
            methods = methods[[nupa]],
            X = X,
            subset = subset,
            cv = cv,
            cf_mat = cf_mat,
            quiet = quiet
          )
        }
      }
      ## Standard stacking
      # Case 1: List of models
    } else {
      if (!inherits(current_model, "ens.learner")) {
        for (i in seq_along(current_model)) {
          sub_element <- current_model[[i]]
          
          if (is.list(sub_element)) {
            subset <- get_subset(sub_element, models, d_mat, z_mat)
            
            nupa_smoother[[i]] <- get_smoother(
              object = sub_element,
              methods = methods[[nupa]],
              X = X,
              subset = subset,
              cv = cv,
              cf_mat = cf_mat,
              quiet = quiet
            )
          }
        }
        # Case 2: Single ensemble model
      } else {
        nupa_smoother <- get_smoother(
          object = current_model,
          methods = methods[[nupa]],
          X = X,
          subset = NULL,
          cv = cv,
          cf_mat = cf_mat,
          quiet = quiet
        )
      }
    }
    results[[nupa]] <- nupa_smoother
  }
  return(results)
}


#' Smoother weights from ensemble learner
#'
#' Extract smoother weights from a short- or standard-stacked ensemble model
#' created by \code{\link[NuisanceParameters]{ensemble_short}} or \code{\link[NuisanceParameters]{ensemble}}.
#'
#' @param object Stacked ensemble learner from \code{\link[NuisanceParameters]{ensemble_short}} or \code{\link[NuisanceParameters]{ensemble}}.
#' @param methods List of methods built via \code{\link[NuisanceParameters]{create_method}}.
#' @param X Covariate matrix of training sample.
#' @param subset Logical vector indicating which observations to use when determining
#'   ensemble weights. If not provided, all observations are used.
#' @param cf_mat Logical matrix with \code{k} columns of indicators representing folds,
#'   typically created by \code{\link[NuisanceParameters]{prep_cf_matrix}}.
#' @param cv Integer specifying the number of cross-validation folds used in estimating
#'   the ensemble model.
#' @param quiet If \code{FALSE}, the method currently computed is printed to the console.
#'
#' @return A matrix of dimension \code{nrow(X)} × \code{nrow(X)} containing ensemble smoother weights.
#'
#' @keywords internal
#'
get_smoother <- function(object,
                         methods,
                         X,
                         subset = NULL,
                         cf_mat, cv,
                         quiet = TRUE) {
  if (is.null(subset)) subset <- rep(TRUE, nrow(X))
  
  ## Short-Stacking
  if (cv == 1) {
    smoother_w <- get_smoother_weights(
      object = object$ens_object, methods = methods, X = X, subset = subset,
      ens_w = object$ens_w, cf_mat = cf_mat, quiet = quiet
    )
    
    ## Standard-Stacking
  } else if (cv > 1) {
    smoother_w <- matrix(0, nrow = nrow(X), ncol = nrow(X))
    
    for (i in seq_len(ncol(cf_mat))) {
      fold <- cf_mat[, i]
      no_cf <- ncol(cf_mat) == 1
      
      X_tr <- X[if (no_cf) subset else (!fold & subset), ]
      X_te <- X[if (no_cf) TRUE else fold, ]
      
      smoother_w_fold <- get_smoother_weights(
        object = object[[i]]$ens_object, methods = methods, 
        ens_w = object[[i]]$ens_w, X = X_tr, Xnew = X_te, quiet = quiet
      )
      
      if (no_cf) {
        smoother_w[, subset] <- smoother_w_fold
      } else {
        smoother_w[fold, !fold & subset] <- smoother_w_fold
      }
    }
  }
  
  return(smoother_w)
}


#' Smoother weights from short-stacked ensemble learner
#'
#' Extract smoother weights from a short-stacked ensemble model created by
#' \code{\link[NuisanceParameters]{ensemble_short}}.
#'
#' @param object Short-stacked ensemble learner from \code{\link[NuisanceParameters]{ensemble_short}}.
#' @param methods List of method models from \code{\link[NuisanceParameters]{create_method}} used
#'  as input for \code{\link[NuisanceParameters]{ensemble_short}}.
#' @param X Covariate matrix of training sample.
#' @param subset Logical vector indicating which observations to use for
#'  determining ensemble weights. If not provided, all observations are used.
#' @param ens_w Ensemble weights to aggregate predictions from different learners.
#'  Must be a vector of length \code{ncol(object$cf_preds)}. Optional.
#' @param cf_mat Logical matrix with \code{k} columns of indicators representing
#'  folds (for example created by \code{\link[NuisanceParameters]{prep_cf_matrix}}).
#' @param quiet If \code{FALSE}, the method currently computed is printed to the console.
#' @param ... Ignore unused arguments.
#'
#' @return A matrix of dimension \code{nrow(X)} × \code{nrow(X)} containing ensemble smoother weights.
#'
#' @method get_smoother_weights ensemble_short
#' @keywords internal
#' @exportS3Method
#'
get_smoother_weights.ensemble_short <- function(object,
                                                methods,
                                                X,
                                                cf_mat,
                                                subset = NULL,
                                                ens_w = NULL,
                                                quiet = TRUE,
                                                ...) {
  
  if (is.null(ens_w)) ens_w <- rep(1 / ncol(object$cf_preds), ncol(object$cf_preds))
  if (is.null(subset)) subset <- rep(TRUE, nrow(X))
  
  ens_models <- object$ens_models
  smoother_w <- matrix(0, nrow = nrow(X), ncol = nrow(X))
  
  ## Smoother extraction
  for (i in seq_len(ncol(cf_mat))) {
    fold <- cf_mat[, i]
    
    # Subset is responsible for T-learner logic
    test_idx <- which(fold)
    train_idx <- which(!fold)
    train_idx_sub <- which(!fold & subset)
    idx_mapping <- match(train_idx_sub, train_idx)
    
    X_tr <- X[train_idx_sub, ]
    X_te <- X[test_idx, ]
    
    # Get T-learner smoother weights (test × train&subset x methods)
    w_array <- get_smoother_weights.ensemble_core(
      object = ens_models[[i]], methods = methods,
      X_tr = X_tr, X_te = X_te, quiet = quiet
    )
    
    # Pad (to test × full train x methods) and insert raw weights
    w_padded <- array(0, dim = c(length(test_idx), length(train_idx), length(methods)))
    w_padded[, idx_mapping, ] <- w_array
    
    smoother_w[fold, !fold] <- agg_array(a = w_padded, w = ens_w)
  }
  
  return(smoother_w)
}


#' Smoother weights from cross-validated ensemble learner
#'
#' Extract ensemble smoother weights for a covariate matrix based on a
#' cross-validated ensemble model created by \code{\link[NuisanceParameters]{ensemble}}.
#'
#' @param object Trained ensemble object from \code{\link[NuisanceParameters]{ensemble}}.
#' @param methods List of methods built via \code{\link[NuisanceParameters]{create_method}}.
#' @param X Covariate matrix of training sample.
#' @param Xnew Covariate matrix of test sample.
#' @param ens_w Ensemble weights to aggregate predictions from different learners.
#'  Must be a vector of length \code{ncol(object$cf_preds)}. Optional.
#' @param quiet If \code{FALSE}, the method currently computed is printed to the console.
#' @param ... Ignore unused arguments.
#'
#' @return A matrix of dimension \code{nrow(Xnew)} × \code{nrow(X)} containing ensemble smoother weights.
#'
#' @method get_smoother_weights ensemble
#' @keywords internal
#' @exportS3Method
#'
get_smoother_weights.ensemble <- function(object,
                                          methods,
                                          X, 
                                          Xnew,
                                          ens_w = NULL,
                                          quiet = TRUE,
                                          ...) {
  w_array <- get_smoother_weights.ensemble_core(
    object = object$ens_models, methods = methods,
    X_tr = X, X_te = Xnew, quiet = quiet
  )
  
  if (length(object$ens_models) > 1) {
    smoother_w <- agg_array(a = w_array, w = ens_w)
  } else {
    smoother_w <- array(w_array[, , 1], dim = dim(w_array)[-3])
  }
  
  return(smoother_w)
}


#' Ensemble smoother weights
#'
#' Extract smoother weights for a (new) covariate matrix based on the fully trained
#' models of an ensemble created by \code{\link[NuisanceParameters]{ensemble_core}}.
#'
#' @param object List of trained models from \code{\link[NuisanceParameters]{ensemble_core}}.
#' @param methods List of methods built via \code{\link[NuisanceParameters]{create_method}}.
#' @param X_tr Covariate matrix of training sample.
#' @param X_te Covariate matrix of test sample.
#' @param quiet If \code{FALSE}, the method currently computed is printed to the console.
#' @param ... Ignore unused arguments.
#'
#' @return A three-dimensional array of dimension
#'   \code{nrow(X_te)} × \code{nrow(X_tr)} × \code{length(methods)} containing smoother weights.
#'   Each slice corresponds to one method, where each row gives the weights
#'   assigned to training outcomes in predicting the respective test observation.
#'
#' @method get_smoother_weights ensemble_core
#' @keywords internal
#' @exportS3Method
#'
get_smoother_weights.ensemble_core <- function(object,
                                               methods,
                                               X_tr, 
                                               X_te,
                                               quiet = FALSE,
                                               ...) {
  w_array <- array(data = 0, dim = c(nrow(X_te), nrow(X_tr), length(object)))
  
  for (i in 1:length(object)) {
    if (isFALSE(quiet)) print(paste0("Smoother weights: ", methods[[i]]$method))

    x_select <- methods[[i]]$x_select
    X_te_sel <- if (is.null(x_select)) X_te else X_te[, x_select]
    
    ## Smoother weights extraction
    w_array[, , i] <- get_smoother_weights(object = object[[i]], Xnew = X_te_sel)
  }
  return(w_array)
}


#' Compute smoother weights from arithmetic mean
#'
#' Generate smoother weights for a test sample using arithmetic mean fitting
#' of the training sample.
#'
#' @param object An object returned by \code{\link[NuisanceParameters]{mean_fit}}.
#' @param Xnew Optional covariate matrix of the test sample.
#'  If not provided, predictions are computed for the training sample.
#' @param ... Ignored arguments.
#'
#' @return A matrix of smoother weights.
#'
#' @keywords internal
#' @method get_smoother_weights mean_fit
#' @exportS3Method
#'
get_smoother_weights.mean_fit <- function(object, Xnew = NULL, ...) {
  if (is.null(Xnew)) {
    n_new <- object$N
  } else {
    n_new <- nrow(Xnew)
  }
  
  S <- matrix(1 / object$N, nrow = n_new, ncol = object$N)
  
  return(S)
}


#' Compute smoother weights from OLS (internal)
#'
#' Generate smoother weights for a test sample from an OLS model.
#'
#' @param object An object returned by \code{\link[NuisanceParameters]{ols_fit}}.
#' @param Xnew Optional covariate matrix of the test sample.
#'  If not provided, predictions are computed for the training sample.
#' @param ... Ignored arguments.
#'
#' @return A matrix of smoother weights.
#'
#' @keywords internal
#' @method get_smoother_weights ols_fit
#' @exportS3Method
#'
get_smoother_weights.ols_fit <- function(object, Xnew = NULL, ...) {
  if (is.null(Xnew)) Xnew <- object$X
  
  X <- add_intercept(object$X)
  Xnew <- add_intercept(Xnew)
  
  # Remove variables dropped due to multi-collinearity
  X <- X[, !is.na(object)]
  Xnew <- Xnew[, !is.na(object)]
  
  S <- Xnew %*% solve(crossprod(X), tol = 2.225074e-308) %*% t(X)
  
  return(S)
}


#' Smoother weights for ridge regression (internal)
#'
#' @description
#' Post-estimation function to extract smoother weights for a ridge regression.
#' The models are estimated via a wrapper of \code{cv.glmnet} in the 
#' \pkg{NuisanceParameters} package with \code{alpha = 0}.
#'
#' @param object An object returned by \code{\link[NuisanceParameters]{ridge_fit}}.
#' @param Xnew Optional covariate matrix of the test sample. If \code{NULL},
#'             the smoother matrix for the training data is returned. If standardization
#'             is used, must be standardized using the same center and scale as \code{X}.
#' @param ... Additional arguments passed to \link{get_smoother_weights}.
#'
#' @return An \eqn{N \times N} smoother matrix.
#'
#' @keywords internal
#' @method get_smoother_weights ridge_fit
#' @exportS3Method
#'
get_smoother_weights.ridge_fit <- function(object, Xnew = NULL, ...) {
  X <- object[["call"]][["x"]]
  Y <- object[["call"]][["y"]]
  
  if (is.null(Xnew)) {
    Xnew <- X
    Xnew <- add_intercept(Xnew)
  } else {
    Xnew <- scale(Xnew, object$x_means, object$x_sds)
    Xnew <- add_intercept(Xnew)
  }
  
  X <- add_intercept(X)
  N <- nrow(X)
  p <- ncol(X) - 1
  
  sd_y <- sqrt(stats::var(Y) * ((N - 1) / N))
  lambda <- (1 / sd_y) * object$lambda.min * N
  
  # Reference: https://stats.stackexchange.com/questions/129179/why-is-glmnet-ridge-regression-giving-me-a-different-answer-than-manual-calculat
  S <- Xnew %*% solve(crossprod(X) + lambda * diag(c(0, rep(1, p)))) %*% t(X)
  
  return(S)
}


#' Smoother weights from k-nearest neighbor algorithm (internal)
#'
#' Extract smoother weights from a k-nearest neighbor algorithm.
#' Each weight equals \code{1 / k} for neighbors and 0 for non-neighbors.
#'
#' @param object Output of \code{\link[NuisanceParameters]{knn_fit}}.
#' @param Xnew Covariate matrix of test sample.
#' If not provided, prediction is done for the training sample.
#' @param ... Ignore unused arguments.
#'
#' @return A matrix of smoother weights.
#'
#' @method get_smoother_weights knn_fit
#' @keywords internal
#' @exportS3Method
#'
get_smoother_weights.knn_fit <- function(object, Xnew = NULL, ...) {
  X <- object$X
  
  if (is.null(Xnew)) Xnew <- X
  k <- if (is.null(object$arguments$k)) 10 else object$arguments$k
  
  # Credit: FastKNN package
  euclidean_dist <- t(apply(Xnew, 1, function(x_i) {
    sqrt(rowSums(sweep(X, 2, x_i)^2))
  }))
  
  get_binary_vec <- function(row) {
    min_indices <- order(row)[1:k]
    binary_vec <- rep(0, length(row))
    binary_vec[min_indices] <- 1
    return(binary_vec)
  }
  
  S <- apply(euclidean_dist, 1, get_binary_vec)
  S <- t(S) / k
  
  return(S)
}


#' Smoother weights for the \code{xgboost_fit} function
#'
#' @description
#' Post-estimation function to extract smoother weights for an XGBoost model
#' estimated via the \code{xgboost_fit} function from the \pkg{NuisanceParameters} package.
#'
#' @param object An object of class \code{xgboost_fit}, i.e., the result of calling \code{xgboost_fit}.
#' @param Xnew Optional covariate matrix of the test sample. If \code{NULL},
#'             the smoother matrix for the training data is returned.
#' @param ... Additional arguments passed to \link{get_smoother_weights}.
#'
#' @return An \eqn{N \times N} smoother matrix.
#'   
#' @method get_smoother_weights xgboost_fit
#' @keywords internal
#' @exportS3Method
get_smoother_weights.xgboost_fit <- function(object,
                                             Xnew = NULL,
                                             ...) {
  X <- object$X
  Y <- object$Y
  
  # Re-route and pass object$booster
  get_smoother_weights.xgb.Booster(
    object = object$booster,
    X = X,
    Y = Y,
    Xnew = Xnew,
    ...
  )
}


## Utils 


#' Determine the Corresponding Data Subset for a Model Element
#'
#' Helper function to identify which column in \code{d_mat} or \code{z_mat} corresponds to
#' a given model element from \code{np_models$Y.hat.d_m} or \code{np_models$Y.hat.z_m}.
#'
#' @param sub_element A model element (e.g., from \code{np_models$Y.hat.d_m} or \code{np_models$Y.hat.z_m}).
#' @param models A list containing ensemble models with components \code{Y.hat.d_m} and \code{Y.hat.z_m}.
#' @param d_mat Logical matrix of treatment indicators.
#' @param z_mat Logical matrix of instrument indicators.
#'
#' @return The corresponding column from \code{d_mat} or \code{z_mat} if a match is found;
#'   otherwise, \code{NULL}.
#'
#' @keywords internal
get_subset <- function(sub_element, models, d_mat, z_mat) {
  # Check for matches in Y.hat.d_m
  idx <- which(vapply(models$Y.hat.d_m, identical, logical(1), sub_element))
  if (length(idx)) {
    return(d_mat[, idx])
  }
  
  # Check for matches in Y.hat.z_m
  idx <- which(vapply(models$Y.hat.z_m, identical, logical(1), sub_element))
  if (length(idx)) {
    return(z_mat[, idx])
  }
  
  NULL
}


#' Aggregate along an array dimension with weights
#'
#' Computes the weighted sum of an array along a specified dimension.
#' By default, aggregates along the second dimension.
#'
#' @param a A numeric array of at least 2 dimensions.
#' @param w A numeric vector of weights. Its length must match the size
#'   of the dimension being aggregated.
#' @param dim Integer. The dimension over which to aggregate
#'   (default is \code{2}).
#'
#' @return A numeric array with the same dimensions as \code{a},
#'   except with the specified dimension removed. For example, if
#'   \code{a} is \eqn{N \times M} and \code{dim = 2}, the result is
#'   a length-\eqn{N} vector. If \code{a} is \eqn{N \times M \times K},
#'   the result is an \eqn{N \times K} matrix.
#'
#' @keywords internal
agg_array <- function(a, w, dim = 2) {
  apply(a, c(1, dim), function(x) sum(x * w))
}


#' Adds an intercept to a matrix
#'
#' \code{\link{add_intercept}} adds an intercept to a matrix.
#'
#' @param mat Any matrix.
#'
#' @return Matrix with intercept.
#'
#' @keywords internal
#'
add_intercept <- function(mat) {
  if (is.null(dim(mat))) {
    mat <- as.matrix(mat, ncol = 1)
  }
  if (all(mat[, 1] == 1)) {
    return(mat)
  } else {
    mat <- cbind(rep(1, nrow(mat)), mat)
    return(mat)
  }
}
