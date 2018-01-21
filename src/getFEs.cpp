#define ARMA_NO_DEBUG
#include <RcppArmadillo.h>


// [[Rcpp::export(name = ".kaczmarz")]]
arma::vec Kaczmarz(const arma::vec &kb,
                   const arma::vec &kw_tilde,
                   const arma::umat &kA,
                   const unsigned int kl,
                   const double kc_alpha_tol,
                   const unsigned int kc_trace) {
  // Auxiliary variables.
  const unsigned int kn = kb.n_rows;
  const unsigned int kK = kA.n_cols;
  
  // Use Kaczmarz to recover the fixed effects.
  arma::vec alpha = arma::zeros(kl);
  unsigned int r = 1;
  while (1) {
    // Check user interrupt.
    Rcpp::checkUserInterrupt();
    
    // Loop through all observations.
    const arma::vec kalpha_old = alpha;
    for (unsigned int i = 0 ; i < kn ; ++i) {
      const arma::uvec kidx = kA.row(i).t();
      arma::vec alpha_tmp = arma::zeros(kl);
      alpha_tmp(kidx).fill((kb(i) - arma::sum(kw_tilde(i) * alpha(kidx))) * 
        kw_tilde(i) / (kK * std::pow(kw_tilde(i), 2)));
      alpha += alpha_tmp;
    }
    
    if (kc_trace > 0) {
      Rcpp::Rcout << "iter= " << r++ 
        << ", crit= " << arma::norm(alpha - kalpha_old) << std::endl;
    }
    
    // Check termination.
    if (arma::norm(alpha - kalpha_old) < kc_alpha_tol) {
      break;
    }
  }
  
  return alpha;
}


// Bugged at the moment.
// // [[Rcpp::export(name = ".closedform")]]
// arma::vec ClosedForm(const arma::vec &kbeta,
//                      const arma::vec &ky,
//                      const arma::mat &kX,
//                      const arma::umat &kD,
//                      const arma::uvec &klvls_k,
//                      const double kc_alpha_tol,
//                      const unsigned int kc_trace) {
//   // Auxiliary variables.
//   const unsigned int kn = ky.n_rows;
//   const unsigned int kK = kD.n_cols;
//   const unsigned int kl = arma::sum(klvls_k);
//   
//   // Sort data to efficiently compute group sums.
//   arma::umat S(kn, kK);
//   for (unsigned int k = 0 ; k < kK ; ++k) {
//     S.col(k) = arma::sort_index(kD.col(k));
//   }
//   
//   // Construct index to subset \alpha.
//   const arma::uvec klvls_k_tmp = arma::cumsum(klvls_k);
//   arma::umat I(kK, 2);
//   I.col(0) = klvls_k_tmp - klvls_k;
//   I.col(1) = klvls_k_tmp - 1;
//   
//   // Compute log(sum(y_k)).
//   arma::vec lsy(kl);
//   for (unsigned int k = 0 ; k < kK ; ++k) {
//     // Sort category k.
//     const arma::vec ky_tmp = ky(S.col(k));
//     arma::uvec d_tmp = kD.col(k);
//     d_tmp = d_tmp(S.col(k));
//     
//     // Sum sorted data.
//     arma::vec sum = arma::zeros(klvls_k(k));
//     unsigned int i = 0;
//     for (unsigned int j = 0 ; j < klvls_k(k) ; ++j) {
//       while (i < kn && d_tmp(i) == j) {
//         sum(j) += ky_tmp(i++);
//       }
//     }
//     
//     // Fill log(sum(y_k)).
//     lsy.subvec(I(k, 0), I(k, 1)) = arma::log(sum);
//   }
//   
//   // Use closed form to recover fixed effects.
//   const arma::vec keta = kX * kbeta;
//   arma::vec alpha = arma::zeros(kl);
//   unsigned int r = 1;
//   while (1) {
//     // Check user interrupt.
//     Rcpp::checkUserInterrupt();
//     
//     // Zig-zaging.
//     const arma::vec kalpha_old = alpha;
//     for (unsigned int k = 0 ; k < kK ; ++k) {
//       // Compute linear predictor - \eta.
//       arma::vec eta = keta;
//       for (unsigned int kk = 0 ; kk < kK ; ++kk) {
//         if (k != kk) {
//           const arma::vec kalpha_k = alpha.subvec(I(kk, 0), I(kk, 1));
//           eta += kalpha_k(kD.col(kk));
//         }
//       } 
//       
//       // Sort category k.
//       const arma::vec keta_tmp = arma::exp(eta(S.col(k)));
//       arma::uvec d_tmp = kD.col(k);
//       d_tmp = d_tmp(S.col(k));
//       
//       // Sum sorted data.
//       arma::vec sum = arma::zeros(klvls_k(k));
//       unsigned int i = 0;
//       for (unsigned int j = 0 ; j < klvls_k(k) ; ++j) {
//         while (i < kn && d_tmp(i) == j) {
//           sum(j) += keta_tmp(i++);
//         }
//       }
//       
//       // Update \alpha.
//       alpha.subvec(I(k, 0), I(k, 1)) = lsy.subvec(I(k, 0), I(k, 1)) - 
//         arma::log(sum);
//     }
//     
//     if (kc_trace > 0) {
//       Rcpp::Rcout << "iter= " << r++ 
//                   << ", crit= " << arma::norm(alpha - kalpha_old) << std::endl;
//     }
//     
//     // Check termination.
//     if (arma::norm(alpha - kalpha_old) < kc_alpha_tol) {
//       break;
//     }
//   }
//   
//   return alpha;
// }