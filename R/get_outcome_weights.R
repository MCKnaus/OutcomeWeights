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
#' @return A list of at least these components:
#' - omega: matrix (number of point estimates x number of estimation units) of outcome weights
#' - treat: the treatment indicator to make it compatible with the cobalt package
#' 
#' @references 
#' Knaus, M. C. (2024). Treatment effect estimators as weighted outcomes, \url{https://arxiv.org/abs/2411.11559}.
#'
#' @export
#'
get_outcome_weights = function(object,...) UseMethod("get_outcome_weights") # makes get_outcome_weights a generic method


#' Outcome weights maker for pseudo-IV estimators.
#'
#' @description This is a generic function taking pseudo-instrument,
#' pseudo-treatment and the transformation matrix as inputs and returning 
#' outcome weights
#'
#' @param Z.tilde Numeric vector of pseudo-instrument outcomes.
#' @param D.tilde Numeric vector of pseudo-treatment.
#' @param T_mat Transformation matrix
#'
#' @return A vector of outcome weights.
#' 
#' @references 
#' Knaus, M. C. (2024). Treatment effect estimators as weighted outcomes, soon on 'arXiv'.
#'
#' @export
#'
pive_weight_maker = function(Z.tilde,D.tilde,T_mat) {
  omega = (t(Z.tilde) %*% D.tilde)^(-1) %*% t(Z.tilde) %*% T_mat
  return(omega)
}




#' \code{summary} method for class `outcome_weights`
#' 
#' Calculates several summary measures of potentially many outcome weights.
#' 
#' @param object \link{get_outcome_weights} object.
#' @param quiet If TRUE, results are passed but not printed.
#' @param digits Number of digits to be displayed. Default 4.
#' @param epsilon Threshold below which in absolute values non-zero but small values should be displayed as < ...
#' @param ... further arguments passed to \code{printCoefmat}
#' 
#' @return 3D-array of dimension 
#' - c("Control","Treated") x
#' - number of point estimates x 
#' - c("Minimum weight","Maximum weight","% Negative","Sum largest 10%","Sum of weights","Sum of absolute weights")
#' 
#' @export
#' 
summary.get_outcome_weights = function(object,
                                       quiet=FALSE,
                                       digits = 4,
                                       epsilon = 0.0001,
                                       ...) {
  
  treat = object$treat
  
  n_pe = nrow(object$omega)
  
  if (isFALSE(quiet) & n_pe > 10) {
    quiet = TRUE
    warning("Weights summary of more then 10 point estimates suppressed to avoid lengthy output.
            Please use the returned array to post-process the results yourself.")
  }
  
  if (is.null(rownames(object$omega))) pe_labels = 1:n_pe
  else pe_labels = rownames(object$omega)
  
  array = array(NA,c(2,n_pe,6),dimnames = list(c("Control","Treated"),pe_labels,
                                            c("Minimum weight","Maximum weight","% Negative","Sum largest 10%","Sum of weights","Sum of absolute weights")))
  
  array[1,,] = summary_weights_rcpp(-object$omega[,treat==0])
  array[2,,] = summary_weights_rcpp(object$omega[,treat==1])
  
  if (isFALSE(quiet)) {
    # Custom function to format and print values with eps notation and footnote
    print_custom_matrix = function(mat, threshold = 0.0001, digits = 4, first_col_width = 25, other_col_width = 10) {
      # Format the matrix and track if "eps" is used
      formatted_matrix = apply(mat, c(1, 2), function(x) {
        if (x == 0) {
          return(sprintf("%.*f", digits, 0))  # Print exact zero as "0.0000"
        } else if (abs(x) < epsilon) {
          if (x > 0) {
            return("epsilon*")   # Small positive
          } else {
            return("-epsilon*")  # Small negative
          }
        } else {
          return(sprintf(paste0("%.", digits, "f"), x))  # Format regular numbers
        }
      })
      
      # Check if "eps" is present in the formatted matrix
      eps_needed = any(formatted_matrix == "epsilon*" | formatted_matrix == "-epsilon*")
      
      # Get row and column names
      row_labels = rownames(mat)
      col_labels = colnames(mat)
      
      # Print column headers with custom widths
      cat(sprintf("%-*s", first_col_width, ""), paste(sprintf("%-*s", other_col_width, col_labels), collapse = " "), "\n")
      
      # Print each row with row labels and aligned column widths
      for (i in 1:nrow(formatted_matrix)) {
        cat(sprintf("%-*s", first_col_width, row_labels[i]), paste(sprintf("%-*s", other_col_width, formatted_matrix[i, ]), collapse = " "), "\n")
      }
      
      # Print footnote if eps values were used
      if (eps_needed) {
        cat("* epsilon < ", toString(epsilon),"\n")
      }
    }
    
    for (e in 1:n_pe) {
      cat("\nWeight summary for estimator", pe_labels[e] ,"\n")
      print_custom_matrix(t(array[, e, ]))
    }
  }
  
  invisible(array)
}
