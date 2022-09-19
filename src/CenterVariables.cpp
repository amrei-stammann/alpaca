#define ARMA_NO_DEBUG
#include <RcppArmadillo.h>

// Method of alternating projections (Halperin)
// [[Rcpp::export(name = "centerVariables")]]
arma::mat CenterVariables(
    const arma::mat& V,
    const arma::vec& w,
    const Rcpp::List& klist,
    const double tol
  ) {
  // Auxiliary variables (fixed)
  const int n = V.n_rows;
  const int K = klist.size();
  const int P = V.n_cols;
  const double sw = arma::accu(w);
  
  // Auxiliary variables (storage)
  double delta, denom, meanj, num, wt;
  int index, iter, j, k, p, t, J, T;
  arma::mat M(n, P);
  arma::vec x(n);
  arma::vec x0(n);
  
  // Halperin projections
  for (p = 0 ; p < P ; ++p) {
    // Center each variable
    x = V.col(p);
    for (iter = 0 ; iter < 100000 ; ++iter) {
      // Check user interrupt
      Rcpp::checkUserInterrupt();
      
      // Store centered vector from the last iteration
      x0 = x;
      
      // Alternate between categories
      for (k = 0 ; k < K ; ++k) {
        // Compute all weighted group means of category 'k' and subtract them
        Rcpp::List jlist = klist[k];
        J = jlist.size();
        for (j = 0 ; j < J ; ++j) {
          // Subset j-th group of category 'k'
          Rcpp::IntegerVector indexes = jlist[j];
          T = indexes.size();
          
          // Compute numerator and denominator of the weighted group mean
          num = 0.0;
          denom = 0.0;
          for (t = 0 ; t < T ; ++t) {
            index = indexes[t];
            wt = w(index);
            num += wt * x(index);
            denom += wt;
          }
          
          // Subtract weighted group mean
          meanj = num / denom;
          for (t = 0 ; t < T ; ++t) {
            index = indexes[t];
            x(index) -= meanj;
          }
        }
      }
      
      // Check convergence
      delta = arma::accu(arma::abs(x - x0) / (1.0 + arma::abs(x0)) % w) / sw;
      if (delta < tol) {
        break;
      }
    }
    M.col(p) = x;
  }
  
  // Return matrix with centered variables
  return M;
}