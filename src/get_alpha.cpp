#define ARMA_NO_DEBUG
#include <RcppArmadillo.h>


// [[Rcpp::export(name = "get.alpha")]]
arma::vec GetAlpha(const arma::vec &kD_alpha,
                   const arma::ivec &klvls_k,
                   const arma::imat &kA,
                   const arma::imat &kB,
                   const double ktol) {
  // Auxiliary variables
  const int kn = kD_alpha.n_rows;
  const int kK = klvls_k.n_rows;
  
  // Start alternating between normal equations
  arma::vec alpha = arma::zeros(arma::accu(klvls_k));
  double crit;
  do {
    // Check user interrupt
    Rcpp::checkUserInterrupt();
    
    // Store \alpha of the previous iteration
    const arma::vec kalpha_old = alpha;
    
    // Solve normal equation 'k'
    int start = 0;
    for (int k = 0 ; k < kK ; ++k) {
      // Last index of \alpha corresponding to category 'k'
      const int kend = start + klvls_k(k) - 1;
      
      // Compute adjusted dependent variable
      int start_k = 0;
      arma::vec b_dots = kD_alpha;
      for (int kk = 0 ; kk < kK ; ++kk) {
        const int kend_k = start_k + klvls_k(kk) - 1;
        if (k != kk) {
          const arma::vec kalpha_k = alpha.subvec(start_k, kend_k);
          for (int l = 0 ; l < kn ; ++l) {
            b_dots(l) -= kalpha_k(kA(l, kk));
          }
        }
        start_k = kend_k + 1;
      }
      
      // Sort category k
      arma::vec v_k(kn);
      arma::ivec a_k(kn);
      for (int i = 0 ; i < kn ; ++i) {
        const int kb = kB(i, k);
        a_k(i) = kA(kb, k);
        v_k(i) = b_dots(kb);
      }
      
      // Group mean of sorted data
      arma::vec alpha_k(klvls_k(k));
      int i = 0;
      for (int j = 0 ; j < klvls_k(k) ; ++j) {
        double sum = 0.0;
        int n = 0;
        while (i < kn && a_k(i) == j) {
          sum += v_k(i);
          ++n;
          ++i;
        }
        alpha_k(j) = sum / n;
      }
      
      // Update \alpha corresponding to category 'k'
      alpha.subvec(start, kend) = alpha_k;
      start = kend + 1;
    }
    
    // Compute termination criterion
    crit = arma::norm(alpha - kalpha_old, 2) / arma::norm(kalpha_old, 2);
  } while (crit >= ktol);
  
  // Return \alpha
  return alpha;
}