#define ARMA_NO_DEBUG
#include <RcppArmadillo.h>


// [[Rcpp::export(name = "pseudo.demeaning")]]
arma::mat PseudoDemeaning(const arma::mat &kV,
                          const arma::vec &kw,
                          const arma::imat &kA,
                          const arma::imat &kB,
                          const arma::ivec &klvls_k,
                          const double ktol) {
  // Auxiliary variables
  const int kn = kV.n_rows;
  const int kp = kV.n_cols;
  const int kk = kA.n_cols;
  
  // Halperin projections (default)
  arma::mat M(kn, kp);
  for (int p = 0 ; p < kp ; ++p) {
    // Pseudo demean variable 'kv'
    arma::vec Mv = kV.col(p);
    double crit;
    do {
      // Check user interrupt
      Rcpp::checkUserInterrupt();
      
      // Alternate between categories
      const arma::vec kMv_old = Mv;
      for (int k = 0 ; k < kk ; ++k) {
        // Sort category k
        arma::vec Mv_k(kn);
        arma::vec w_k(kn);
        arma::ivec a_k(kn);
        for (int i = 0 ; i < kn ; ++i) {
          const int kb = kB(i, k);
          a_k(i) = kA(kb, k);
          Mv_k(i) = Mv(kb);
          w_k(i) = kw(kb);
        }
        
        // Pseudo-demean sorted data
        arma::vec fac_k(klvls_k(k));
        int i = 0;
        for (int j = 0 ; j < klvls_k(k) ; ++j) {
          double num = 0.0;
          double denom = 0.0;
          while (i < kn && a_k(i) == j) {
            num += w_k(i) * Mv_k(i);
            denom += std::pow(w_k(i), 2);
            ++i;
          }
          fac_k(j) = num / denom;
        }
        
        // Sort back and substract
        for (int i = 0 ; i < kn ; ++i) {
          Mv(i) -= kw(i) * fac_k(kA(i, k));
        }
      }
      
      // Compute termination criterion
      crit = arma::norm(Mv - kMv_old, 2) / arma::norm(kMv_old, 2);
    } while (crit >= ktol);
    
    // Add column to M
    M.col(p) = Mv;
  }
  
  // Return matrix with demeaned variable
  return M;
}