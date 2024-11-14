#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;

// [[Rcpp::export]]
arma::mat scaled_Ztildex_maker(arma::sp_mat& alpha, const arma::colvec& res_z, const arma::colvec& res_d)  {
  
  int n_train = alpha.n_cols; // get the dimension of training data
  int n_test = alpha.n_rows; // get the dimensions of test data
  
  arma::mat Ztildex(n_train,n_test); // set up an empty matrix
  arma::mat scaled_Ztildex(n_train,n_test); // set up an empty matrix 
  
  alpha = trans(alpha);

  for(int i = 0; i < n_test; i++) {
    arma::colvec alphai(alpha.col(i));
    arma::colvec demeaned_Zres = res_z - arma::dot(res_z, alphai); // this is R'M in the theory formula (demeaned residual)
    Ztildex.col(i) = demeaned_Zres % alphai; // this is R'M diag(alpha)
    // Calculate the sum of the element-wise product of Ztildex's i-th column and Dres, i.e. the thing in ()^-1 in theory formula
    double denominator = arma::sum(Ztildex.col(i) % res_d);
    // Divide the i-th column of Ztildex by this sum and assign it to scaled_Ztildex's i-th column
    scaled_Ztildex.col(i) = Ztildex.col(i) / denominator;
    // Now multiply the respective column with the rows of trans_subtracted to get the first column of the weight matrix
    Rcpp::checkUserInterrupt();
  }
  
  // Perform matrix multiplication between 'scaled_Ztildex' and the result of the subtraction
  scaled_Ztildex = trans(scaled_Ztildex);
  return scaled_Ztildex;
}
