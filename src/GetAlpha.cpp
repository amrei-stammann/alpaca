#define ARMA_NO_DEBUG
#include <RcppArmadillo.h>


// [[Rcpp::export(name = "getAlpha")]]
arma::field<arma::vec> GetAlpha(const arma::vec& pi,
                                const Rcpp::List& klist,
                                const double tol) {
  // Auxiliary variables (fixed)
  const int n = pi.n_rows;
  const int K = klist.size();
  
  // Auxiliary variables (storage)
  double crit, denom, num, sum;
  int iter, j, k, kk, t, J, T;
  arma::vec y(n);
  
  // Generate starting guess
  arma::field<arma::vec> Alpha(K);
  for (k = 0 ; k < K ; ++k) {
    Rcpp::List jlist = klist[k];
    J = jlist.size();
    Alpha(k) = arma::zeros(J);
  }
  
  // Start alternating between normal equations
  arma::field<arma::vec> Alpha0(arma::size(Alpha));
  for (iter = 0 ; iter < 10000 ; ++iter) {
    // Check user interrupt
    Rcpp::checkUserInterrupt();
    
    // Store \alpha_{0} of the previous iteration
    Alpha0 = Alpha;
    
    // Solve normal equations of category k
    for (k = 0 ; k < K ; ++k) {
      // Compute adjusted dependent variable
      y = pi;
      for (kk = 0 ; kk < K ; ++kk) {
        if (kk != k) {
          Rcpp::List jlist = klist[kk];
          J = jlist.size();
          for (j = 0 ; j < J ; ++j) {
            Rcpp::IntegerVector indexes = jlist[j];
            T = indexes.size();
            for (t = 0 ; t < T ; ++t) {
              y(indexes[t]) -= Alpha(kk)(j);
            }
          }
        }
      }
      
      // Compute group mean
      Rcpp::List jlist = klist[k];
      J = jlist.size();
      arma::vec alpha(J);
      for (j = 0 ; j < J ; ++j) {
        // Subset the j-th group of category k
        Rcpp::IntegerVector indexes = jlist[j];
        T = indexes.size();
        
        // Compute group sum
        sum = 0.0;
        for (t = 0 ; t < T ; ++t) {
          sum += y(indexes[t]);
        }
        
        // Store group mean
        alpha(j) = sum / T;
      }
      
      
      // Update \alpha_{k}
      Alpha(k) = alpha;
    }
    
    // Compute termination criterion and check convergence
    num = 0.0;
    denom = 0.0;
    for (k = 0 ; k < K ; ++k) {
      num += arma::accu(arma::pow(Alpha(k) - Alpha0(k), 2));
      denom += arma::accu(arma::pow(Alpha0(k), 2));
    }
    crit = std::sqrt(num / denom);
    if (crit < tol) {
      break;
    }
  }
  
  // Return \alpha
  return Alpha;
}