#' Smoother weights method
#'
#' @description This is a generic method for getting smoother weights.
#' It calculates the smoother weights for objects created by other packages.
#' See get_smoother_weight.<compatible_class> in the package documentation for compatible functions.
#'
#' @param object An object obtained from other packages.
#' @param ... Additional arguments specific to object class implementations.
#'            See the documentation for which objects require which additional arguments.
#'
#' @return An \eqn{N \times N} smoother matrix.
#'
#' @export
#'
get_smoother_weights <- function(object, ...) UseMethod("get_smoother_weights")


#' Smoother weights for the \code{estimatr::lm_robust} function
#'
#' @description
#' Post-estimation function to extract smoother weights for a linear regression
#' estimated via the \code{lm_robust} function from the \pkg{estimatr} package.
#'
#' @param object An object of class \code{lm_robust}, i.e., the result of calling \code{estimatr::lm_robust()}.
#' @param X Covariate matrix of the training sample.
#' @param Xnew Optional covariate matrix of the test sample. If \code{NULL},
#'             the smoother matrix for the training data is returned.
#' @param ... Additional arguments passed to \link{get_smoother_weights}.
#'
#' @return An \eqn{N \times N} smoother matrix.
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#' n <- 100
#' X <- matrix(rnorm(n * 3), ncol = 3)
#' Y <- 2 * X[, 1] + 3 * X[, 2] - 1 * X[, 3] + rnorm(n)
#'
#' # Fit linear model
#' object <- estimatr::lm_robust(Y ~ 0 + X)
#'
#' # Get fitted values
#' preds <- object$fitted.values
#'
#' # Get smoother matrix
#' S <- get_smoother_weights(object, X = X, Xnew = X)
#'
#' # Check equivalence
#' all.equal(as.numeric(S %*% Y), as.numeric(preds))
#' }
#'
#' @method get_smoother_weights lm_robust
#' @export
#'
get_smoother_weights.lm_robust <- function(object, X, Xnew, ...) {
  if (is.null(Xnew)) Xnew <- X
  has_intercept <- "(Intercept)" %in% names(object$coefficients)

  if (has_intercept) {
    X <- add_intercept(X)
    Xnew <- add_intercept(Xnew)
  }

  S <- Xnew %*% solve(crossprod(X), tol = 2.225074e-308) %*% t(X)

  return(S)
}


#' Smoother weights for the \code{stats::lm} function
#'
#' @description
#' Post-estimation function to extract smoother weights for a linear regression
#' estimated via the \code{lm} function from the \pkg{stats} package.
#'
#' @param object An object of class \code{lm}, i.e., the result of calling \code{stats::lm()}.
#' @param X Covariate matrix of the training sample.
#' @param Xnew Optional covariate matrix of the test sample. If \code{NULL},
#'             the smoother matrix for the training data is returned.
#' @param ... Additional arguments passed to \link{get_smoother_weights}.
#'
#' @return An \eqn{N \times N} smoother matrix.
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#' n <- 100
#' X <- matrix(rnorm(n * 3), ncol = 3)
#' Y <- 2 * X[, 1] + 3 * X[, 2] - 1 * X[, 3] + rnorm(n)
#'
#' # Fit linear model
#' object <- lm(Y ~ 0 + X)
#'
#' # Get fitted values
#' preds <- object$fitted.values
#'
#' # Get smoother matrix
#' S <- get_smoother_weights(object, X = X, Xnew = X)
#'
#' # Check equivalence
#' all.equal(as.numeric(S %*% Y), as.numeric(preds))
#' }
#'
#' @method get_smoother_weights lm
#' @export
#'
get_smoother_weights.lm <- function(object, X, Xnew, ...) {
  if (is.null(Xnew)) Xnew <- X
  has_intercept <- "(Intercept)" %in% names(object$coefficients)

  if (has_intercept) {
    X <- add_intercept(X)
    Xnew <- add_intercept(Xnew)
  }

  S <- Xnew %*% solve(crossprod(X), tol = 2.225074e-308) %*% t(X)

  return(S)
}


#' Smoother weights for the \code{plasso::cv.plasso} function
#'
#' @description
#' Post-estimation function to extract smoother weights for a post-Lasso regression
#' estimated via the \code{cv.plasso} function from the \pkg{plasso} package.
#'
#' @param object An object of class \code{cv.plasso}, i.e., the result of 
#'               calling \code{plasso::cv.plasso()}.
#' @param Xnew Optional covariate matrix of the test sample. If \code{NULL},
#'             the smoother matrix for the training data is returned.
#' @param ... Additional arguments passed to \link{get_smoother_weights}.
#'
#' @return An \eqn{N \times N} smoother matrix.
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#' n <- 100
#' X <- matrix(rnorm(n * 3), ncol = 3)
#' Y <- 2 * X[, 1] + 3 * X[, 2] - 1 * X[, 3] + rnorm(n)
#'
#' # Fit post-lasso model
#' object <- plasso::cv.plasso(x = X, y = Y)
#'
#' # Get predictions
#' preds <- predict(object, newx = X, type = "response", s = "optimal", se_rule = 0)$plasso
#'
#' # Get smoother matrix
#' S <- get_smoother_weights(object, Xnew = X)
#'
#' # Check equivalence
#' all.equal(as.numeric(S %*% Y), as.numeric(preds))
#' }
#'
#' @references
#' Knaus, M. C. (2024). Treatment effect estimators as weighted outcomes, \url{https://arxiv.org/abs/2411.11559}.
#'
#' @method get_smoother_weights cv.plasso
#' @export
#'
get_smoother_weights.cv.plasso <- function(object, Xnew = NULL, ...) {
  if (is.null(Xnew)) Xnew <- object$x

  X <- add_intercept(object$x)
  Xnew <- add_intercept(Xnew)

  colnames(X)[1] <- "(Intercept)"
  colnames(Xnew) <- colnames(X)

  # Handles no intercept
  Xact <- X[, object$names_pl, drop = FALSE]
  Xactnew <- Xnew[, object$names_pl, drop = FALSE]

  S <- Xactnew %*% solve(crossprod(Xact), tol = 2.225074e-308) %*% t(Xact)

  return(S)
}


#' Smoother weights for the \code{hdm::rlasso} function
#'
#' @description
#' Post-estimation function to extract smoother weights for a post-Lasso regression
#' estimated via the \code{rlasso} function from the \pkg{hdm} package.
#'
#' @param object An object of class \code{rlasso}, i.e., the result of 
#'               calling \code{hdm::rlasso()}.
#' @param Xnew Optional covariate matrix of the test sample. If \code{NULL},
#'             the smoother matrix for the training data is returned.
#' @param ... Additional arguments passed to \link{get_smoother_weights}.
#'
#' @return An \eqn{N \times N} smoother matrix.
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#' n <- 100
#' X <- matrix(rnorm(n * 3), ncol = 3)
#' Y <- 2 * X[, 1] + 3 * X[, 2] - 1 * X[, 3] + rnorm(n)
#'
#' # Fit post-lasso model
#' object <- hdm::rlasso(x = X, y = Y, post = TRUE)
#'
#' # Get predictions
#' preds <- predict(object, newdata = X)
#'
#' # Get smoother matrix
#' S <- get_smoother_weights(object, Xnew = X)
#'
#' # Check equivalence
#' all.equal(as.numeric(S %*% Y), as.numeric(preds))
#' }
#' 
#' @references 
#' Knaus, M. C. (2024). Treatment effect estimators as weighted outcomes, \url{https://arxiv.org/abs/2411.11559}.
#'
#' @method get_smoother_weights rlasso
#' @export
#'
get_smoother_weights.rlasso <- function(object, Xnew = NULL, ...) {
  if (!object$options$post) stop("Smoother matrix requires post-lasso specification.")
  
  X <- object$model
  if (is.null(Xnew)) Xnew <- X
  active_vars <- which(object$index)

  if ("(Intercept)" %in% names(object$coefficients)) {
    X <- add_intercept(X)
    Xnew <- add_intercept(Xnew)
    colnames(X)[1] <- "(Intercept)"
    colnames(Xnew) <- colnames(X)
    active_vars <- c("(Intercept)", names(object$beta)[active_vars])
  }

  Xact <- X[, active_vars, drop = FALSE]
  Xactnew <- Xnew[, active_vars, drop = FALSE]

  S <- Xactnew %*% solve(crossprod(Xact), tol = 2.225074e-308) %*% t(Xact)

  return(S)
}


#' Smoother weights for the \code{grf::regression_forest} function
#'
#' @description
#' Post-estimation function to extract smoother weights for a random forest model
#' estimated via the \code{regression_forest} function from the \pkg{grf} package.
#'
#' @param object An object of class \code{regression_forest}, i.e.,
#'               the result of calling \code{grf::regression_forest()}.
#' @param Xnew Optional covariate matrix of the test sample. If \code{NULL},
#'             the smoother matrix for the training data is returned.
#' @param ... Additional arguments passed to \link{get_smoother_weights}.
#'
#' @return An \eqn{N \times N} smoother matrix.
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#' n <- 100
#' X <- matrix(rnorm(n * 3), ncol = 3)
#' Y <- 2 * X[, 1] + 3 * X[, 2] - 1 * X[, 3] + rnorm(n)
#'
#' # Fit regression forest
#' object <- grf::regression_forest(X = X, Y = Y)
#'
#' # Get predictions
#' preds <- predict(object, newdata = X)$predictions
#'
#' # Get smoother matrix
#' S <- get_smoother_weights(object)
#'
#' # Check equivalence
#' all.equal(as.numeric(S %*% Y), as.numeric(preds))
#' }
#'
#' @references 
#' Athey, S., Tibshirani, J., & Wager, S. (2019). Generalized random forest. The Annals of Statistics, 47(2), 1148-1178.
#' 
#' Knaus, M. C. (2024). Treatment effect estimators as weighted outcomes, \url{https://arxiv.org/abs/2411.11559}.
#' 
#' @method get_smoother_weights regression_forest
#' @export
#'
get_smoother_weights.regression_forest <- function(object, Xnew = NULL, ...) {
  if (is.null(Xnew)) Xnew <- object$X.orig

  S <- grf::get_forest_weights(object, newdata = Xnew)
  S <- as.matrix(S)

  return(S)
}


#' Smoother weights for the \code{drf::drf} function
#'
#' @description
#' Post-estimation function to extract smoother weights for a distributional random
#' forest model estimated via the \code{drf} function from the \pkg{drf} package.
#'
#' @param object An object of class \code{drf}, i.e., the result of calling \code{drf::drf()}.
#' @param Xnew Optional covariate matrix of the test sample. If \code{NULL},
#'             the smoother matrix for the training data is returned.
#' @param ... Additional arguments passed to \link{get_smoother_weights}.
#'
#' @return An \eqn{N \times N} smoother matrix.
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#' n <- 100
#' X <- matrix(rnorm(n * 3), ncol = 3)
#' Y <- 2 * X[, 1] + 3 * X[, 2] - 1 * X[, 3] + rnorm(n)
#'
#' # Fit distributional forest
#' object <- drf::drf(X = X, Y = Y)
#'
#' # Get predictions
#' preds <- predict(object, newdata = X, functional = "mean")[["mean"]]
#'
#' # Get smoother matrix
#' S <- get_smoother_weights(object)
#'
#' # Check equivalence
#' all.equal(as.numeric(S %*% Y), as.numeric(preds))
#' }
#'
#' @method get_smoother_weights drf
#' @export
#'
get_smoother_weights.drf <- function(object, Xnew = NULL, ...) {
  if (is.null(Xnew)) Xnew <- object$X.orig

  S <- as.matrix(predict(object, newdata = Xnew)$weights)

  return(S)
}


#' Smoother weights for the \code{xgboost::xgb.train} function
#'
#' @description
#' Post-estimation function to extract smoother weights for an XGBoost model
#' estimated via the \code{xgb.train} function from the \pkg{xgboost} package.
#'
#' @param object An object of class \code{xgb.Booster}, i.e., the result of 
#'               calling \code{xgboost::xgb.train()}.
#' @param X Covariate matrix of the training sample.
#' @param Y Vector of outcomes of the training sample.
#' @param Xnew Optional covariate matrix of the test sample. If \code{NULL},
#'             the smoother matrix for the training data is returned.
#' @param ... Additional arguments passed to \link{get_smoother_weights}.
#'
#' @return An \eqn{N \times N} smoother matrix.
#'
#' @note For smoother extraction, the XGBoost model must be trained with:
#' \code{alpha = 0}, \code{subsample = 1}, \code{max_delta_step = 0}.
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#' n <- 100
#' X <- matrix(rnorm(n * 3), ncol = 3)
#' Y <- 2 * X[, 1] + 3 * X[, 2] - 1 * X[, 3] + rnorm(n)
#'
#' # Set the parameters
#' nrounds <- 100
#' params <- list(alpha = 0, subsample = 1, max_delta_step = 0)
#'
#' # Prepare the data
#' d_train <- xgboost::xgb.DMatrix(data = as.matrix(X), label = Y)
#' d_test <- xgboost::xgb.DMatrix(data = as.matrix(X))
#'
#' # Fit XGBoost
#' object <- xgboost::xgb.train(data = d_train, nrounds = nrounds, params = params)
#'
#' # Get predictions
#' preds <- predict(object, newdata = d_test)
#'
#' # Get smoother matrix
#' S <- get_smoother_weights(object, X, Y)
#'
#' # Check equivalence
#' all.equal(as.numeric(S %*% Y), as.numeric(preds))
#' }
#'
#' @method get_smoother_weights xgb.Booster
#' @export
#'
get_smoother_weights.xgb.Booster <- function(object,
                                             X = NULL,
                                             Y = NULL,
                                             Xnew = NULL,
                                             ...) {
  if (is.null(Xnew)) Xnew <- X
  params <- attr(object, "params")

  # Params for smoother to be extracted
  supported_vals <- list(
    alpha          = 0,
    subsample      = 1,
    max_delta_step = 0
  )

  supported_names <- names(supported_vals)
  check_difference <- function(k) {!identical(params[[k]], supported_vals[[k]])}
  
  # Check against supported values
  diff <- supported_names[vapply(supported_names, check_difference, logical(1))]
  
  if (length(diff) > 0) {
    stop(
      paste0(
        "Smoothers are available for a subset of hyperparameters.\n",
        "These parameter(s) differ:  ",
        paste(diff, collapse = ", "), "\n",
        "Type ?get_smoother_weights.xgb.Booster to see the supported values."
      ),
      call. = FALSE
    )
  }

  dtrain <- xgboost::xgb.DMatrix(data = as.matrix(X), label = Y)
  dtest <- xgboost::xgb.DMatrix(data = as.matrix(Xnew))

  S <- get_xgboost_weights(model = object, dtrain = dtrain, dtest = dtest)
  S <- as.matrix(S$S_test)

  return(S)
}


#' Smoother weights for the \code{ranger::ranger} function
#'
#' @description
#' Post-estimation function to extract smoother weights for a random forest model
#' estimated via the \code{ranger} function from the \pkg{ranger} package.
#'
#' @param object An object of class \code{ranger}, i.e., the result of 
#'               calling \code{ranger::ranger()}.
#' @param X Optional covariate matrix of the training sample; 
#'          used if not stored in the fitted object.
#' @param Xnew Covariate matrix of test sample.
#'             If not provided, prediction is done for the training sample.
#' @param ... Additional arguments passed to \link{get_smoother_weights}.
#'
#' @return An \eqn{N \times N} smoother matrix.
#' 
#' @examples
#' \dontrun{
#' set.seed(123)
#' n <- 100
#' X <- matrix(rnorm(n * 3), ncol = 3, dimnames = list(NULL, 1:3))
#' Y <- 2 * X[, 1] + 3 * X[, 2] - 1 * X[, 3] + rnorm(n)
#'
#' # Fit random forest
#' object <- ranger::ranger(x = X, y = Y, keep.inbag = TRUE)
#'
#' # Get predictions
#' preds <- predict(object, data = X)$predictions
#'
#' # Get smoother matrix
#' S <- get_smoother_weights(object, X = X)
#'
#' # Check equivalence
#' all.equal(as.numeric(S %*% Y), as.numeric(preds))
#' }
#'
#' @method get_smoother_weights ranger
#' @export
#'
get_smoother_weights.ranger <- function(object, 
                                        X = NULL, 
                                        Xnew = NULL, 
                                        ...) {
  # Credit: With some adjustments, this function is borrowed from 
  # https://github.com/maxi-tb22/DualML
  
  if (is.null(object$inbag.counts)) {
    stop("Set 'keep.inbag=TRUE' in ranger::ranger()")
  }
  
  X <- if (!is.null(X)) X else object$call$x
  if (is.null(X)) stop("Please provide the covariate matrix X used for training.")
  if (is.null(Xnew)) Xnew <- X
  
  leaves_train <- predict(object, X, type = "terminalNodes", predict.all = TRUE)$predictions
  leaves_test  <- predict(object, Xnew,  type = "terminalNodes", predict.all = TRUE)$predictions
  
  n_trees <- ncol(leaves_train)
  n_train <- nrow(X)
  n_test  <- nrow(Xnew)
  
  S <- matrix(0, nrow = n_test, ncol = n_train)
  
  for (t in seq_len(n_trees)) {
    inbag_counts_t <- object$inbag.counts[[t]]
    leaf_train_t   <- leaves_train[, t]
    leaf_test_t    <- leaves_test[, t]
    
    for (j in seq_len(n_test)) {
      lt <- leaf_test_t[j]
      idx <- which(leaf_train_t == lt & inbag_counts_t > 0)
      if (length(idx)) {
        counts <- inbag_counts_t[idx]
        denom  <- sum(counts)
        # Add this tree's contribution (will normalize by n_trees at end)
        S[j, idx] <- S[j, idx] + counts / denom
      }
    }
  }
  # Average over trees
  S <- S / n_trees
  return(S)
}
