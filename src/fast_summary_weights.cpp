#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;

// [[Rcpp::export]]
arma::mat summary_weights_rcpp(arma::mat weights) {
  
  int n = weights.n_rows;
  int nw = weights.n_cols;
  
  arma::mat summary(n, 6);
  arma::rowvec temp(nw);
  int p90 = static_cast<int>(0.9 * nw);
  
  for (int i = 0; i < n; i++) {
    //Rcpp::Rcout << "Row " << i << std::endl;
    temp = arma::sort(weights.row(i));
    summary(i, 0) = temp.min();
    summary(i, 1) = temp.max();
    
    // Share with negative weights
    if (summary(i, 0) < 0) {
      for (int j = 0; j < nw; j++) {
        if (temp(j) >= 0) {
          summary(i, 2) = static_cast<double>(j) / nw;
          break;
        }
      }
    } else {
      summary(i, 2) = 0.0;
    }
    
    // Sum of largest 10%
    summary(i, 3) = arma::sum(temp.subvec(p90, nw - 1));
    
    // Sum of weights
    summary(i, 4) = arma::sum(temp);
    
    // Sum of absolute weights
    summary(i, 5) = arma::sum(arma::abs(temp));
    Rcpp::checkUserInterrupt();
  }
  
  return summary;
}
