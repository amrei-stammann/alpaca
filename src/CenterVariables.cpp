#define ARMA_NO_DEBUG
#include <RcppArmadillo.h>

// Method of alternating projections (Halperin)
// [[Rcpp::export(name = "centerVariables")]]
arma::mat CenterVariables(const arma::mat& V,
                          const arma::vec& w,
                          const Rcpp::List& klist,
                          const double tol) {
  // Auxiliary variables (fixed)
  const int n = V.n_rows;
  const int K = klist.size();
  const int P = V.n_cols;
  
  // Auxiliary variables (storage)
  double adjtol, c, crit, delta, delta0, denom, meanj, num, wt;
  int index, iter, j, k, p, t, J, T;
  arma::mat M = V;
  arma::vec x0(n);
  
  // Halperin projections
  for (p = 0 ; p < P ; ++p) {
    // Center each variable
    delta0 = 2.0 * arma::norm(M.col(p), 2);
    for (iter = 0 ; iter < 10000 ; ++iter) {
      // Check user interrupt
      Rcpp::checkUserInterrupt();
      
      // Store vector from the last iteration
      x0 = M.col(p);
      
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
            num += wt * M(index, p);
            denom += std::pow(wt, 2);
          }
          
          // Subtract weighted group mean
          meanj = num / denom;
          for (t = 0 ; t < T ; ++t) {
            index = indexes[t];
            M(index, p) -= w(index) * meanj;
          }
        }
      }
      
      // Compute termination criteria and check convergence
      delta = arma::norm(M.col(p) - x0, 2);
      c = delta / delta0;
      adjtol = std::sqrt(1.0 + arma::accu(arma::pow(x0, 2))) * tol;
      crit = delta / (1.0 - c);
      if (crit <= adjtol) {
        break;
      }
      
      // Store \delta from the previous iteration
      delta0 = delta;
    }
  }
  
  // Return matrix with centered variables
  return M;
}