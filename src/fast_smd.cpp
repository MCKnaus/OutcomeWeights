#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;

// [[Rcpp::export]]
arma::cube smd_rcpp(arma::mat x, arma::rowvec treat, arma::mat weights, arma::mat targets) {
  
  int n = x.n_rows;
  int d = x.n_cols;
  int nw = weights.n_rows;
  
  arma::cube results(d, 5, nw);
  arma::mat variances(2, d);
  
  arma::rowvec temp_w(n);
  arma::rowvec temp_t(d);
  arma::mat temp(5, d);
  
  for(int i = 0; i < 2; i++) {
    arma::mat temp_x = x.rows(arma::find(treat == i));
    variances.row(i) = arma::var(temp_x);
  }
  
  arma::rowvec nd_sd = arma::sqrt(arma::mean(variances));
  
  for(int i = 0; i < nw; i++) {
    // Rcpp::Rcout << "Row " << i << std::endl;
    temp_w = weights.row(i);
    temp_t = targets.row(i);
    temp.row(0) = (temp_w % (treat - 1)) * x;
    temp.row(1) = (temp_w % treat) * x;
    temp.row(2) = (temp.row(1) - temp.row(0)) / nd_sd;
    temp.row(3) = (temp.row(0) - temp_t) / nd_sd;
    temp.row(4) = (temp.row(1) - temp_t) / nd_sd;
    results.slice(i) = temp.t();
    
    Rcpp::checkUserInterrupt();
  }
  
  return results;
}

bool greater20(double i) { return (i > 20); }
bool greater10(double i) { return (i > 10); }
bool greater5(double i) { return (i > 5); }

// [[Rcpp::export]]
arma::cube summary_smd_rcpp(arma::cube cbr) {
  
  int n = cbr.n_slices;
  double nx = cbr.n_rows; // double because we need fractional shares
  int ncol = cbr.n_cols;
  
  arma::cube summary(6, 2, n);
  arma::mat temp_s(nx, ncol);
  arma::vec temp_b(nx);
  arma::vec temp_t(nx * 2);
  arma::mat temp(6, 2);
  arma::uvec target_id = {3, 4};
  
  for(int i = 0; i < n; i++) {
    temp_s = arma::abs(cbr.slice(i));
    temp_b = temp_s.col(2);
    temp_t = arma::vectorise(temp_s.cols(target_id));
    temp(0, 0) = arma::max(temp_s.col(2));
    temp(0, 1) = arma::max(temp_t);
    temp(1, 0) = arma::mean(temp_s.col(2));
    temp(1, 1) = arma::mean(temp_t);
    temp(2, 0) = arma::median(temp_s.col(2));
    temp(2, 1) = arma::median(temp_t);
    temp(3, 0) = std::count_if(temp_b.begin(), temp_b.end(), greater20) / nx;
    temp(3, 1) = std::count_if(temp_t.begin(), temp_t.end(), greater20) / (nx * 2);
    temp(4, 0) = std::count_if(temp_b.begin(), temp_b.end(), greater10) / nx;
    temp(4, 1) = std::count_if(temp_t.begin(), temp_t.end(), greater10) / (nx * 2);
    temp(5, 0) = std::count_if(temp_b.begin(), temp_b.end(), greater5) / nx;
    temp(5, 1) = std::count_if(temp_t.begin(), temp_t.end(), greater5) / (nx * 2);
    summary.slice(i) = temp;
    Rcpp::checkUserInterrupt();
  }
  
  return summary;
}
