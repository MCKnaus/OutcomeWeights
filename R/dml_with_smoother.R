### This is a package internal implementation of Double ML.


#' This function runs Double ML estimators with outcome smoothers.
#'
#' Existing Double ML implementations are too general to easily extract smoother matrices
#' required to be compatible with the get_forest_weights() method. This motivates yet 
#' another Double ML implementation.
#'
#' @param Y Numeric vector containing the outcome variable.
#' @param D Optional binary treatment variable.
#' @param X Covariate matrix with N rows and p columns.
#' @param Z Optional binary instrumental variable.
#' @param estimators String (vector) indicating which estimators should be run.
#' Current menu: c("PLR","PLR_IV","AIPW_ATE","Wald_AIPW")
#' @param smoother Indicate which smoother to be used for nuisance parameter estimation.
#' Currently only available option \code{"honest_forest"} from the \pkg{grf} package.
#' @param n_cf_folds Number of cross-fitting folds. Default is 5.
#' @param n_reps Number of repetitions of cross-fitting. Default is 1.
#' @param ... Options to be passed to smoothers.
#' 
#' @return A list with three entries:
#' - \code{results}: a list storing the results, influence functions, and score functions of each estimator
#' - \code{NuPa.hat}: a list storing the estimated nuisance parameters and the outcome smoother matrices
#' 
#' @export
#' 
dml_with_smoother = function(Y,D,X,Z=NULL,
                             estimators = c("PLR","PLR_IV","AIPW_ATE","Wald_AIPW"),
                             smoother = "honest_forest", 
                             n_cf_folds = 5,
                             n_reps=1,
                             ...) {
  # Sanity checks
  supported_estimators = c("PLR","PLR_IV","AIPW_ATE","Wald_AIPW")
  not_supported = estimators[!estimators %in% supported_estimators]
  if (length(not_supported) > 0) {
    stop(paste("Error: The following specified estimators are not supported:", 
               paste(not_supported, collapse = ", ")))}
  supported_smoothers = c("honest_forest")
  not_supported = smoother[!smoother %in% supported_smoothers]
  if (length(not_supported) > 0) {
    stop(paste("Error: The following specified smoothers are not supported:", 
               paste(not_supported, collapse = ", ")))}
  if (any(c("PLR_IV", "Wald_AIPW") %in% estimators) && is.null(Z)) {
    stop("Error: Z cannot be NULL when using either 'PLR_IV' or 'Wald_AIPW' as an estimator.")
  }
  
  # Extract useful information
  N = length(Y)
  n_estimators = length(estimators)
  
  # Find required NuPas
  NuPa = character(0)
  if ("PLR" %in% estimators) NuPa = c(NuPa, "Y.hat", "D.hat")
  if ("PLR_IV" %in% estimators) NuPa = c(NuPa, "Y.hat", "D.hat", "Z.hat")
  if ("AIPW_ATE" %in% estimators) NuPa = c(NuPa, "Y.hat.d", "D.hat")
  if ("Wald_AIPW" %in% estimators) NuPa = c(NuPa, "Y.hat.z", "D.hat.z", "Z.hat")
  NuPa = unique(NuPa)
  
  # Estimate required nuisance parameters
  if (smoother == "honest_forest") {
  NuPa.hat = NuPa_honest_forest(NuPa = NuPa,
                             X,Y,D,Z,
                             n_cf_folds=n_cf_folds,
                             n_reps=n_reps,
                             ...) }
  
  # Intialize empty
  dml_plr = dml_PLR_IV = dml_AIPW_ATE = dml_Wald_AIPW = 
    "This estimator was not run."
  
  # Run the specified DML estimators
  if ("PLR" %in% estimators) {
    D.hat = NuPa.hat$predictions$D.hat
    Y.hat = NuPa.hat$predictions$Y.hat
    psi.a = -(D - D.hat)^2
    psi.b = (Y - Y.hat) * (D - D.hat)
    dml_plr = dml_inference(psi.a,psi.b)
  }
  if ("PLR_IV" %in% estimators) {
    D.hat = NuPa.hat$predictions$D.hat
    Y.hat = NuPa.hat$predictions$Y.hat
    Z.hat = NuPa.hat$predictions$Z.hat
    psi.a = -(D - D.hat) * (Z - Z.hat)
    psi.b = (Y - Y.hat) * (Z - Z.hat) 
    dml_PLR_IV = dml_inference(psi.a,psi.b)
  }
  if ("AIPW_ATE" %in% estimators) {
    D.hat = NuPa.hat$predictions$D.hat
    Y.hat.d0 = NuPa.hat$predictions$Y.hat.d0
    Y.hat.d1 = NuPa.hat$predictions$Y.hat.d1
    psi.a = matrix(-1,N,n_reps)
    psi.b = Y.hat.d1 - Y.hat.d0 + D * (Y - Y.hat.d1) / D.hat - (1 - D) * (Y - Y.hat.d0) / (1-D.hat)
    dml_AIPW_ATE = dml_inference(psi.a,psi.b)
  }
  if ("Wald_AIPW" %in% estimators) {
    Z.hat = NuPa.hat$predictions$Z.hat
    D.hat.z0 = NuPa.hat$predictions$D.hat.z0
    D.hat.z1 = NuPa.hat$predictions$D.hat.z1
    Y.hat.z0 = NuPa.hat$predictions$Y.hat.z0
    Y.hat.z1 = NuPa.hat$predictions$Y.hat.z1
    psi.a = -( D.hat.z1 - D.hat.z0 + Z * (D - D.hat.z1) / Z.hat - (1 - Z) * (D - D.hat.z0) / (1-Z.hat) )
    psi.b = Y.hat.z1 - Y.hat.z0 + Z * (Y - Y.hat.z1) / Z.hat - (1 - Z) * (Y - Y.hat.z0) / (1-Z.hat)
    dml_Wald_AIPW = dml_inference(psi.a,psi.b)
  }
  
  list_results = list(
    "PLR" = dml_plr,
    "PLR_IV" = dml_PLR_IV,
    "AIPW_ATE" = dml_AIPW_ATE,
    "Wald_AIPW" = dml_Wald_AIPW )
  
  list_data = list(
    "Y" = Y,
    "D" = D,
    "Z" = Z,
    "X" = X )
  
  list_nums = list(
    "N" = N,
    "n_estimators" = n_estimators,
    "n_reps" = n_reps,
    "n_cf_folds" = n_cf_folds
  )
  
  output = list("results" = list_results,
                "NuPa.hat" = NuPa.hat,
                "data" = list_data,
                "numbers" = list_nums)
  
  class(output) = c("dml_with_smoother")
  
  return(output)
}



#' Outcome weights for the \code{\link[grf]{causal_forest}} function
#'
#' @description Post-estimation command to extract outcome weights for double ML
#' run with an outcome smoother.
#'
#' @param object An object of class \code{dml_with_smoother}, i.e. the result of running \code{\link{dml_with_smoother}}.
#' @param ... Pass potentially generic \link{get_outcome_weights} options.
#' @param all_reps If \code{TRUE}, outcomes weights of each repetitions passed. Default \code{TRUE}.
#'
#' @return 
#' - If \code{all_reps == TRUE}: Matrix of dimension (number of estimators x N) containing the outcome weights.
#' - If \code{all_reps == FALSE}: list with two components:
#'      - omega: Matrix of dimension (number of estimators x N) containing the outcome weights.
#'      - omega_all_reps: A list containing the outcome weights of each repetition.
#'
#' @export
#' 
get_outcome_weights.dml_with_smoother = function(object,..., all_reps = FALSE) {
  
  # Preps
  N = object$numbers$N
  n_reps = object$numbers$n_reps
  omega = estimator_names = NULL
  if (n_reps == 1) all_reps = FALSE
  
  # Go through those estimators called and extract their weights
  if (!is.character(object$results$PLR)) {
    Z.tilde = D.tilde = object$data$D - object$NuPa.hat$predictions$D.hat
    omega_plr = NULL
    for (r in 1:n_reps) {
      T_mat = diag(N) - object$NuPa.hat$smoothers$S[r,,]
      omega_plr = rbind(omega_plr, pive_weight_maker(Z.tilde[,r], D.tilde[,r], T_mat) )
    }
    if(isFALSE(all.equal(as.numeric(omega_plr %*% object$data$Y),
                         as.numeric(object$results$PLR$TaPa[,1])))){
      warning("Estimated PLR using weights differ from original estimates.") }
    estimator_names = c(estimator_names, "PLR") 
    omega = rbind(omega,colMeans(omega_plr))  }
  
  if (!is.character(object$results$PLR_IV)){
    Z.tilde = object$data$Z - object$NuPa.hat$predictions$Z.hat
    D.tilde = object$data$D - object$NuPa.hat$predictions$D.hat
    omega_plriv = NULL
    for (r in 1:n_reps) {
      T_mat = diag(N) - object$NuPa.hat$smoothers$S[r,,]
      omega_plriv = rbind(omega_plriv, pive_weight_maker(Z.tilde[,r], D.tilde[,r], T_mat) )
    }
    if(isFALSE(all.equal(as.numeric(omega_plriv %*% object$data$Y),
                         as.numeric(object$results$PLR_IV$TaPa[,1])))){
      warning("Estimated PLR-IV using weights differ from original estimates.") }
    estimator_names = c(estimator_names, "PLR-IV") 
    omega = rbind(omega,colMeans(omega_plriv))  }
  
  if (!is.character(object$results$AIPW_ATE)){
    Z.tilde = D.tilde = matrix(1,N,n_reps)
    lambda1 = object$data$D / object$NuPa.hat$predictions$D.hat
    lambda0 = (1 - object$data$D) / (1 - object$NuPa.hat$predictions$D.hat)
    omega_aipw = NULL
    for (r in 1:n_reps) {
      T_mat = object$NuPa.hat$smoothers$S.d1[r,,] - object$NuPa.hat$smoothers$S.d0[r,,] + 
        lambda1[,r] * (diag(N) - object$NuPa.hat$smoothers$S.d1[r,,]) - 
        lambda0[,r] * (diag(N) - object$NuPa.hat$smoothers$S.d0[r,,])
      omega_aipw = rbind(omega_aipw, pive_weight_maker(Z.tilde[,r], D.tilde[,r], T_mat) )
    }
    if(isFALSE(all.equal(as.numeric(omega_aipw %*% object$data$Y),
                         as.numeric(object$results$AIPW_ATE$TaPa[,1])))){
      warning("Estimated AIPW-ATE using weights differ from original estimates.") }
    estimator_names = c(estimator_names, "AIPW-ATE") 
    omega = rbind(omega,colMeans(omega_aipw))  }
  
  if (!is.character(object$results$Wald_AIPW)){
    lambdaz1 = object$data$Z / object$NuPa.hat$predictions$Z.hat
    lambdaz0 = (1 - object$data$Z) / (1 - object$NuPa.hat$predictions$Z.hat)
    Z.tilde = matrix(1,N,n_reps)
    D.tilde = object$NuPa.hat$predictions$D.hat.z1 - object$NuPa.hat$predictions$D.hat.z0 + 
      lambdaz1 * (object$data$D - object$NuPa.hat$predictions$D.hat.z1) - 
      lambdaz0 * (object$data$D - object$NuPa.hat$predictions$D.hat.z0)
    omega_waipw = NULL
    for (r in 1:n_reps) {
      T_mat = object$NuPa.hat$smoothers$S.z1[r,,] - object$NuPa.hat$smoothers$S.z0[r,,] + 
        lambdaz1[,r] * (diag(N) - object$NuPa.hat$smoothers$S.z1[r,,]) - 
        lambdaz0[,r] * (diag(N) - object$NuPa.hat$smoothers$S.z0[r,,])
      omega_waipw = rbind(omega_waipw, pive_weight_maker(Z.tilde[,r], D.tilde[,r], T_mat) )
    }
    if(isFALSE(all.equal(as.numeric(omega_waipw %*% object$data$Y),
                         as.numeric(object$results$Wald_AIPW$TaPa[,1])))){
      warning("Estimated Wald-AIPW using weights differ from original estimates.") }
    estimator_names = c(estimator_names, "Wald-AIPW") 
    omega = rbind(omega,colMeans(omega_waipw))  }

  # Pass results
  if (isFALSE(all_reps)) {rownames(omega) = estimator_names; return(omega)}
  else {
    list_all_weights = list(
      "PLR" = omega_plr,
      "PLR_IV" = omega_plriv,
      "AIPW_ATE" = omega_aipw,
      "Wald_AIPW" = omega_waipw
    )
    output = list(
      "omega" = omega,
      "omega_all_reps" = list_all_weights
    )
    return(output)
  }
}



#' Results and influence functions for linear double ML estimators.
#'
#' @param psi.a psi.a component of linear Neyman orthogonal score
#' @param psi.b psi.a component of linear Neyman orthogonal score
#'
#' @return List of three components:
#' - \code{TaPa} Matrix storing estimate, standard error, t-value and p-value of target parameter of each repetition in the rows
#' - \code{IF} Matrix storing the influence function of the target parameter of each repetition in the columns
#' - \code{score} Array of dimension n_reps x N x 3 where the matrix in the first dimension stores the score in the first, 
#' psi.a in the second, and psi.b in the third column
#'
#' @keywords internal
#'
dml_inference = function(psi.a,psi.b) {
  # Extract important information
  N = nrow(psi.a)
  n_reps = ncol(psi.a)
  
  # Initialize output containers
  TaPa = matrix(NA,n_reps,4)
  colnames(TaPa) = c("Estimate","SE","t","p")
  IF = matrix(NA,N,n_reps)
  score = array(NA,c(n_reps,N,3))
  dimnames(score) = list( NULL, NULL, third_dim_names = c("psi","psi.a","psi.b") )

  # Loop over replications
  for (r in 1:n_reps) {
    theta = -sum(psi.b[,r]) / sum(psi.a[,r])
    psi = theta * psi.a[,r] + psi.b[,r]
    IF[,r] = - psi / mean(psi.a[,r])
    se = sqrt(var(IF[,r]) / N)
    t = theta / se
    p = 2 * pt(abs(t),N,lower.tail = FALSE)
    TaPa[r,] = c(theta,se,t,p)
    score[r,,] = cbind(psi,psi.a[,r],psi.b[,r])
  }
  
  list("TaPa"=TaPa, "IF"=IF, "score"=score)
}



#' \code{summary} method for class \code{\link{dml_with_smoother}}
#'
#' @param object Object of class \code{\link{dml_with_smoother}}.
#' @param contrast Tests the differences between the coefficients.
#' @param quiet If TRUE, results are passed but not printed.
#' @param ... further arguments passed to \code{printCoefmat}
#'
#' @export
#'
summary.dml_with_smoother = function(object,
                                     contrast=FALSE,
                                     quiet=FALSE,
                                     ...) {
  # Little function to repeatedly create results
  results_maker = function(theta, IF) {
    N = length(IF)
    se = sqrt(var(IF) / N)
    t = theta / se
    p = 2 * pt(abs(t),N,lower.tail = FALSE)
    return(c(theta,se,t,p))
  }
  
  # Extract results for each active estimator
  results = IF = estimator_names = NULL
  if (!is.character(object$results$PLR)) {
    Psi = rowMeans(object$results$PLR$IF)
    results = rbind(results, results_maker(mean(object$results$PLR$TaPa[,1]),Psi))
    IF = cbind(IF, Psi)
    estimator_names = c(estimator_names, "PLR")}
  if (!is.character(object$results$PLR_IV)){
    Psi = rowMeans(object$results$PLR_IV$IF)
    results = rbind(results, results_maker(mean(object$results$PLR_IV$TaPa[,1]),Psi))
    IF = cbind(IF, Psi)
    estimator_names = c(estimator_names, "PLR-IV")}
  if (!is.character(object$results$AIPW_ATE)){
    Psi = rowMeans(object$results$AIPW_ATE$IF)
    results = rbind(results, results_maker(mean(object$results$AIPW_ATE$TaPa[,1]),Psi))
    IF = cbind(IF, Psi)
    estimator_names = c(estimator_names, "AIPW-ATE")}
  if (!is.character(object$results$Wald_AIPW)){
    Psi = rowMeans(object$results$Wald_AIPW$IF)
    results = rbind(results, results_maker(mean(object$results$Wald_AIPW$TaPa[,1]),Psi))
    IF = cbind(IF, Psi)
    estimator_names = c(estimator_names, "Wald-AIPW")}
  
  # Name and print results
  colnames(results) = c("Estimate","SE","t","p")
  rownames(results) = estimator_names
  if(isFALSE(quiet) & isFALSE(contrast)) printCoefmat(results,has.Pvalue = TRUE)
  
  
  if(isTRUE(contrast)) {
    # Number of rows in the matrix
    rows = nrow(results)
    
    # Prepare a vector to store the differences
    differences = SEs = numeric()
    labels = character()
    
    # Nested loops to calculate differences
    for (i in 1:(rows-1)) {
      for (j in (i+1):rows) {
        diff_value = results[j, 1] - results[i, 1]  # Difference of estimates
        diff_IF = IF[,j] - IF[,i]
        diff_label = paste(rownames(results)[j], "-", rownames(results)[i])
        
        # Store the results
        differences = c(differences, diff_value)
        SEs = c(SEs,sqrt(var(diff_IF)/length(diff_IF)))
        labels = c(labels, diff_label)
      }
    }
    
    t = differences / SEs
    p = 2 * pt(abs(t),length(diff_IF),lower.tail = FALSE)
    
    results = cbind(differences,SEs,t,p)
    
    # Assign names to the differences
    rownames(results) = labels
    colnames(results) = c("Estimate","SE","t","p")
    
    if(isFALSE(quiet)) printCoefmat(results,has.Pvalue = TRUE)
  } 
  
  invisible(results)
}



#' \code{plot} method for class \code{\link{dml_with_smoother}}
#'
#' @param x Object of class \code{\link{dml_with_smoother}}.
#' @param ... Pass generic \code{\link[base]{plot}} options.
#' @param alpha Significance level for confidence intervals (default 0.05).
#' @param contrast Shows the differences between the coefficients.
#' 
#' @import ggplot2
#'
#' @export
#'
#' @method plot dml_with_smoother
plot.dml_with_smoother = function(x,..., 
                                  alpha=0.05,
                                  contrast = FALSE) {
  # This is only here to avoid the "no visible global function definition"
  Estimate = Estimator = Lower = Upper = NULL
  
  # Extract results
  results = summary(x,quiet=TRUE,contrast = contrast)
  
  # Z value for 95% confidence
  N = length(x$data$Y)
  t_critical = qt(1 - alpha / 2, N-1)
  
  # Calculate Confidence Intervals
  results = as.data.frame(results)
  results$Estimator = rownames(results)
  results$Estimator = factor(results$Estimator, levels = rownames(results))
  results$Lower = results[,1] - t_critical * results[,2]
  results$Upper = results[,1] + t_critical * results[,2]
  
  # Creating the plot with dodge points and error lines without caps
  p = ggplot(results, aes(x = Estimator, y = Estimate)) +
    geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.1, color = "black") +  # Plain error lines without caps
    geom_point(size = 4, color = "blue") +  # Dodge points to avoid overlap
    geom_hline(yintercept = 0, linetype = "solid", color = "black") +  # Horizontal line at y = 0
    labs(x = "Estimator",
         y = "Estimate") +
    theme_minimal() +
    theme(  # Additional theming to adjust text and layout
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 12),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10)
    )
   if (isTRUE(contrast)) p = p + coord_flip() + 
                              labs(x = "Estimators compared",
                                   y = "Difference")
  
  # Print the plot
  print(p)
}


