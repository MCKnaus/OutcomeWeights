#' Compute Smoother Matrix for a Single Tree
#' 
#' @description Computes the smoother matrix \eqn{S} for a single decision tree.
#'
#' @param train_leaf_ids_tree Numeric vector. Leaf indices for the training data in the tree.
#' @param test_leaf_ids_tree Numeric vector. Leaf indices for the test data in the tree.
#' @param lambda Numeric. The L2 regularization term (lambda) used in the XGBoost model.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{S_train}: The smoother matrix for the training set (n_train x n_train).
#'   \item \code{S_test}: The smoother matrix for the test set (n_test x n_train).
#' }
#' @keywords internal
#' 
create_S_from_tree <- function(train_leaf_ids_tree,
                               test_leaf_ids_tree,
                               lambda) {
  n_train <- length(train_leaf_ids_tree)
  n_test <- length(test_leaf_ids_tree)
  
  S_train <- matrix(0, nrow = n_train, ncol = n_train)
  S_test <- matrix(0, nrow = n_test, ncol = n_train)
  
  leaf_counts <- table(train_leaf_ids_tree)
  unique_leaves <- as.numeric(names(leaf_counts))
  
  for (leaf in unique_leaves) {
    count <- as.numeric(leaf_counts[as.character(leaf)])
    val <- 1 / (count + lambda)
    
    train_idx <- which(train_leaf_ids_tree == leaf)
    test_idx <- which(test_leaf_ids_tree == leaf)
    
    S_train[train_idx, train_idx] <- val
    S_test[test_idx, train_idx] <- val
  }
  
  return(list(S_train = S_train, S_test = S_test))
}


#' Create Smoother Matrix for a Single Boosted Tree
#' 
#' @description Computes the smoother matrix contributions for a specific boosted tree step,
#' applying corrections based on the previous ensemble state (\code{prev_S_train}).
#'
#' @param train_leaf_ids_tree Numeric vector. Leaf indices for the training data in the tree.
#' @param test_leaf_ids_tree Numeric vector. Leaf indices for the test data in the tree.
#' @param prev_S_train Matrix or NULL. The accumulated smoother matrix from previous trees.
#'                     If NULL, this is treated as the first tree.
#' @param lambda Numeric. The L2 regularization parameter.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{S_train}: The smoother matrix for the training set (n_train x n_train).
#'   \item \code{S_test}: The smoother matrix for the test set (n_test x n_train).
#' }
#' @keywords internal
#' 
create_S_from_single_boosted_tree <- function(train_leaf_ids_tree,
                                              test_leaf_ids_tree,
                                              prev_S_train,
                                              lambda) {
  S <- create_S_from_tree(
    train_leaf_ids_tree, test_leaf_ids_tree, lambda
  )
  
  if (is.null(prev_S_train)) {
    # first tree: just normal tree
    return(S)
  }
  
  leaf_sums <- rowsum(prev_S_train, group = train_leaf_ids_tree)
  leaf_counts <- as.numeric(table(train_leaf_ids_tree))
  leaf_corrections <- leaf_sums / (leaf_counts + lambda)
  
  unique_nodes <- as.numeric(rownames(leaf_sums))
  
  train_map <- match(train_leaf_ids_tree, unique_nodes)
  test_map  <- match(test_leaf_ids_tree, unique_nodes)
  
  S_train_correction <- leaf_corrections[train_map, , drop = FALSE]
  S_test_correction  <- leaf_corrections[test_map, , drop = FALSE]
  
  S_train <- S$S_train - S_train_correction
  S_test  <- S$S_test  - S_test_correction
  
  return(list(S_train = S_train, S_test = S_test))
}


#' Iterative Smoother Matrix Construction from XGBoost Regressor
#' 
#' @description Iterates through all trees in an XGBoost model to 
#' construct the final smoother matrix.
#'
#' @param model An object of class \code{xgb.Booster}.
#' @param leaf_indices_train Matrix. Leaf indices for the training data (rows=samples, cols=trees).
#' @param leaf_indices_test Matrix. Leaf indices for the test data.
#' @param base_score Numeric or character. The model’s initial prediction (0 or "mean").
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{S_train}: The smoother matrix for the training set (n_train x n_train).
#'   \item \code{S_test}: The smoother matrix for the test set (n_test x n_train).
#' }
#' @keywords internal
#' 
create_S_from_gbtregressor <- function(model,
                                       leaf_indices_train,
                                       leaf_indices_test,
                                       base_score) {
  n_train <- nrow(leaf_indices_train)
  n_test <- nrow(leaf_indices_test)
  
  # RR: Check for aliases and apply a default [xgboost specificity]
  get_param <- function(params, aliases, default) {
    for (alias in aliases) {
      if (!is.null(params[[alias]])) { return(params[[alias]]) }}
    return(default)
  }
  
  # Fetch eta and lambda from model
  params <- attr(model, "params")
  eta <- get_param(params, c("eta", "learning_rate"), 0.3)
  lambda <- get_param(params, c("lambda", "reg_lambda"), 1)
  
  if (base_score == 0) {
    S_train_curr <- matrix(0, nrow = n_train, ncol = n_train)
    S_test_curr <- matrix(0, nrow = n_test, ncol = n_train)
  } else if (base_score == "mean") {
    S_train_curr <- matrix(1 / n_train, nrow = n_train, ncol = n_train)
    S_test_curr <- matrix(1 / n_train, nrow = n_test, ncol = n_train)
  }
  
  number_of_trees <- ncol(leaf_indices_train)
  
  for (col in 1:number_of_trees) {
    current_tree_indices_train <- leaf_indices_train[, col]
    current_tree_indices_test <- leaf_indices_test[, col]
    
    if (col == 1 && base_score == 0) {
      S <- create_S_from_single_boosted_tree(
        current_tree_indices_train,
        current_tree_indices_test,
        NULL,
        lambda
      )
    } else {
      S <- create_S_from_single_boosted_tree(
        current_tree_indices_train,
        current_tree_indices_test,
        S_train_curr,
        lambda
      )
    }
    
    S_train_curr <- S_train_curr + eta * S$S_train
    S_test_curr <- S_test_curr + eta * S$S_test

  }
  
  return(list(S_train = S_train_curr, S_test = S_test_curr))
}


#' Extract XGBoost Smoother Weights
#' 
#' @description Main function to extract the smoother matrix from a trained 
#' XGBoost model. It validates model parameters and computes 
#' the smoother matrices for training and test sets.
#'
#' @param model An object of class \code{xgb.Booster}. The trained model.
#' @param dtrain An \code{xgb.DMatrix} object. The training data used to build the model.
#' @param dtest An \code{xgb.DMatrix} object. The test data.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{S_train}: The smoother matrix for the training set (n_train x n_train).
#'   \item \code{S_test}: The smoother matrix for the test set (n_test x n_train).
#' }
#' @keywords internal
#' 
get_xgboost_weights <- function(model,
                                dtrain,
                                dtest) {
  params <- attr(model, "params")
  
  label_mean <- mean(xgboost::getinfo(dtrain, "label"))
  base_score <- NULL

  if (!is.null(params$base_score)) {
    if (isTRUE(all.equal(params$base_score, label_mean))) {
      base_score <- "mean"
    } else if (identical(params$base_score, 0)) {
      base_score <- 0
    }
  }
  
  if (is.null(base_score)) stop("base_score must be 0 or mean(Y) for xgboost call; DoubleML call requires 0.")
  
  leaf_indices_train <- predict(model, dtrain, predleaf = TRUE)
  leaf_indices_test <- predict(model, dtest, predleaf = TRUE)
  
  smoothers <- create_S_from_gbtregressor(
    model, leaf_indices_train, leaf_indices_test, base_score
  )
  
  return(list(S_train = smoothers$S_train, S_test = smoothers$S_test))
}
