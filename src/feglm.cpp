#define ARMA_NO_DEBUG
#include <RcppArmadillo.h>


/*
 * Distributions
 */


double lbeta(const double ka,
             const double kb) {
  return std::lgamma(ka) + std::lgamma(kb) - std::lgamma(ka + kb);
}


double ldbinom(const double kk,
               const double kn,
               const double kp) {
  return (kn - kk) * std::log(1.0 - kp) + kk * std::log(kp) - 
    lbeta(kk + 1.0, kn - kk + 1.0) - std::log(kn + 1.0);
}


double ldpois(const double kk,
              const double klambda) {
  return - klambda + kk * std::log(klambda) - std::lgamma(kk + 1.0);
}


/*
 * GLM's
 */


arma::vec InverseLink(const arma::vec &keta,
                      const unsigned int kfamily) {
  // Auxiliary variables.
  const double keps = std::numeric_limits<double>::epsilon() * 10.0;
  
  // Compute inverse link of \eta.
  const unsigned int kn = keta.n_rows;
  arma::vec mu(kn);
  if (kfamily == 0) {
    // --- Logit ---
    mu = 1.0 / (1.0 + arma::exp(- keta));
    
    // Safeguard mu.
    mu(arma::find(mu < keps)).fill(keps);
    mu(arma::find(mu > 1.0 - keps)).fill(1.0 - keps);
  } else if (kfamily == 2) {
    // --- Poisson ---
    mu = arma::exp(keta);
    
    // Safeguard mu.
    mu(arma::find(mu < keps)).fill(keps);
  }
  
  return mu;
}


arma::vec DetaDmu(const arma::vec &kmu,
                  const unsigned int kfamily) {
  // Compute derivative of link function wrt \mu.
  const unsigned int kn = kmu.n_rows;
  arma::vec detadmu(kn);
  if (kfamily == 0) {
    // --- Logit ---
    detadmu = 1.0 / (kmu % (1.0 - kmu));
  } else if (kfamily == 2) {
    // --- Poisson ---
    detadmu = 1.0 / kmu;
  }
  
  return detadmu;
}


arma::vec Variance(const arma::vec &kmu,
                   const unsigned int kfamily) {
  // Compute variance.
  const unsigned int kn = kmu.n_rows;
  arma::vec v(kn);
  if (kfamily == 0) {
    // --- Logit ---
    v = kmu % (1.0 - kmu);
  } else if(kfamily == 2) {
    // --- Poisson ---
    v = kmu;
  }
  
  return v;
}


void Transforming(const arma::vec &kbeta,
                  const arma::vec &kD_alpha,
                  const arma::vec &ky,
                  const arma::mat &kX,
                  const unsigned int kfamily,
                  arma::vec &y_tilde,
                  arma::mat &X_tilde,
                  arma::vec &w_tilde) {
  // Compute linear predictor - \eta.
  const arma::vec keta = kX * kbeta + kD_alpha;

  // Compute \mu.
  const arma::vec kmu = InverseLink(keta, kfamily);
  
  // Compute \frac{\partial \eta}{\partial \mu}.
  const arma::vec kdeta_dmu = DetaDmu(kmu, kfamily); 
  
  // Compute V(\mu).
  const arma::vec kv = Variance(kmu, kfamily);
  
  // Compute \tilde{w}.
  w_tilde = arma::sqrt(1.0 / (arma::pow(kdeta_dmu, 2) % kv));
  
  // Compute \tilde{y}.
  y_tilde = w_tilde % (ky - kmu) % kdeta_dmu;
  
  // Compute \tilde{X}.
  X_tilde = kX.each_col() % w_tilde;
}


arma::vec PseudoDemeaning(const arma::vec &kv,
                          const arma::vec &kw_tilde,
                          const arma::umat &kD,
                          const arma::umat &kS,
                          const arma::uvec &klvls_k,
                          const double kc_pseudo_tol) {
  // Auxiliary variables.
  const unsigned int kn = kv.n_rows;
  const unsigned int kK = kD.n_cols;
  
  // Apply pseudo-demeaning.
  arma::vec Mv = kv;
  while (1) {
    // Check user interrupt.
    Rcpp::checkUserInterrupt();
    
    // Start pseudo-demeaning.
    const arma::vec kMv_old = Mv;
    for (unsigned int k = 0 ; k < kK ; ++k) {
      // Sort category k.
      const arma::vec kMv_tmp = Mv(kS.col(k));
      arma::uvec d_tmp = kD.col(k);
      d_tmp = d_tmp(kS.col(k));
      const arma::vec kw_tilde_tmp = kw_tilde(kS.col(k));
      
      // Pseudo-demean sorted data.
      arma::vec nom = arma::zeros(klvls_k(k));
      arma::vec denom = arma::zeros(klvls_k(k));
      unsigned int i = 0;
      for (unsigned int j = 0 ; j < klvls_k(k) ; ++j) {
        while (i < kn && d_tmp(i) == j) {
          nom(j) += kw_tilde_tmp(i) * kMv_tmp(i);
          denom(j) += std::pow(kw_tilde_tmp(i++), 2);
        }
      }
      
      // Sort back and substract.
      const arma::vec kfac = nom / denom;
      Mv -= kw_tilde % kfac(kD.col(k));
    }
    
    // Check termination.
    if (arma::norm(Mv - kMv_old) < kc_pseudo_tol) {
      break;
    }
  }
  
  return Mv;
}


double LogLikelihood(const arma::vec &kbeta,
                     const arma::vec &kD_alpha,
                     const arma::vec &ky,
                     const arma::mat &kX,
                     const unsigned int kfamily) {
  // Auxiliary variables.
  const unsigned int kn = ky.n_rows;
  const double kmax = std::numeric_limits<double>::max();
  
  // Compute linear predictor - \eta.
  const arma::vec keta = kX * kbeta + kD_alpha;
  
  // Compute \mu.
  const arma::vec kmu = InverseLink(keta, kfamily);
  
  // Compute sum of the loglikelihood.
  double L = 0.0;
  if (kfamily == 0) {
    // --- Logit ---
    for (unsigned int i = 0 ; i < kn ; ++i) {
      const double kx = ldbinom(ky(i), 1.0, kmu(i));
      L += std::isfinite(kx) ? kx : - kmax;
    }
  } else if (kfamily == 2) {
    // --- Poisson ---
    for (unsigned int i = 0 ; i < kn ; ++i) {
      const double kx = ldpois(ky(i), kmu(i));
      L += std::isfinite(kx) ? kx : - kmax;
    }
  }
  
  return L;
}


double xlogy(const double kx,
             const double ky) {
  return kx > 0.0 ? kx * std::log(ky) : 0.0;
}


double Deviance(const arma::vec &kbeta,
                const arma::vec &kD_alpha,
                const arma::vec &ky,
                const arma::mat &kX,
                const unsigned int kfamily) {
  // Auxiliary variables.
  const unsigned int kn = ky.n_rows;
  const double keps = std::numeric_limits<double>::epsilon();
  const double kmax = std::numeric_limits<double>::max();
  
  // Compute linear predictor - \eta.
  const arma::vec keta = kX * kbeta + kD_alpha;
  
  // Compute \mu.
  const arma::vec kmu = InverseLink(keta, kfamily);
  
  // Compute deviance.
 double d = 0.0;
  if (kfamily == 0) {
    for (unsigned int i = 0 ; i < kn ; ++i) {
      if (ky(i) == 1) {
        d -= 2.0 * std::log(kmu(i));
      } else {
        const double kx = kmu(i) == 1.0 ? kmu(i) - keps : kmu(i);
        d -= 2.0 * std::log(1.0 - kx);
      }
    }
  } else if (kfamily == 2) {
    for (unsigned int i = 0 ; i < kn ; ++i) {
      const double kx = 2.0 * (xlogy(ky(i), ky(i) / kmu(i)) - ky(i) + kmu(i));
      d += std::isfinite(kx) ? kx : kmax;
    }
  }
  
  return d;
}


// [[Rcpp::export(name = ".feglm")]]
Rcpp::List FEGLM(const arma::vec &ky,
                 const arma::mat &kX,
                 const arma::umat &kD,
                 const arma::uvec &klvls_k,
                 const arma::vec &kbeta_start,
                 const arma::vec &kD_alpha_start,
                 const unsigned int kfamily,
                 const Rcpp::List &kctrl) {
  // Extract parameters related to control.
  const double kc_grad_tol = Rcpp::as<double>(kctrl["grad.tol"]);
  const double kc_step_tol = Rcpp::as<double>(kctrl["step.tol"]);
  const double kc_dev_tol = Rcpp::as<double>(kctrl["dev.tol"]);
  const double kc_pseudo_tol = Rcpp::as<double>(kctrl["pseudo.tol"]);
  const double kc_rho_tol = Rcpp::as<double>(kctrl["rho.tol"]);
  const unsigned int kc_iter_max = Rcpp::as<unsigned int>(kctrl["iter.max"]);
  const unsigned int kc_trace = Rcpp::as<unsigned int>(kctrl["trace"]);
  
  // Auxiliary variables.
  const unsigned int kK = kD.n_cols;
  const unsigned int kn = ky.n_rows;
  const unsigned int kp = kX.n_cols;
  arma::umat S(kn, kK);
  for (unsigned int k = 0 ; k < kK ; ++k) {
    S.col(k) = arma::sort_index(kD.col(k));
  }
  
  // Compute initial deviance.
  arma::vec beta = kbeta_start;
  arma::vec D_alpha = kD_alpha_start;
  double dev = Deviance(beta, D_alpha, ky, kX, kfamily);
  double ll = LogLikelihood(beta, D_alpha, ky, kX, kfamily);
  if (kc_trace > 1) {
    Rcpp::Rcout << "------------------------------" << std::endl;
    Rcpp::Rcout << "Initial values:" << std::endl;
    Rcpp::Rcout << "n= " << kn << std::endl;
    Rcpp::Rcout << "Deviance= " << dev << std::endl;
    Rcpp::Rcout << "ll= " << ll << std::endl;
  }
  
  // Start IWLS.
  arma::vec b(kn);
  arma::vec g(kp);
  arma::vec w_tilde(kn);
  arma::mat G(kn, kp);
  arma::mat H(kp, kp);
  unsigned int conv;
  unsigned int iter = 1;
  while (1) {
    // Check user interrupt.
    Rcpp::checkUserInterrupt();
    
    // Trace.
    if (kc_trace > 0) {
      Rcpp::Rcout << "--- Iteration= " << iter << " ---" << std::endl;
      if (kc_trace > 1) {
        Rcpp::Rcout << "1. Transform Variables" << std::endl;
      }
    }
    
    // Compute transformed variables.
    arma::vec y_tilde(kn);
    arma::mat X_tilde(kn, kp);
    Transforming(beta, D_alpha, ky, kX, kfamily, y_tilde, X_tilde, w_tilde);
    
    // Trace.
    if (kc_trace > 1) {
      Rcpp::Rcout << "2. Pseudo-Demean Variables" << std::endl;
    }
    
    // Compute pseudo-demeaned variables.
    const arma::vec kMy_tilde = PseudoDemeaning(y_tilde, w_tilde, kD, S,
                                                klvls_k, kc_pseudo_tol);
    arma::mat MX_tilde(kn, kp);
    for (unsigned int j = 0 ; j < kp ; ++j) {
      MX_tilde.col(j) = PseudoDemeaning(X_tilde.col(j), w_tilde, kD, S, klvls_k,
                   kc_pseudo_tol);
    }
    
    // Compute update of the structural parameters.
    g = MX_tilde.t() * kMy_tilde;
    H = MX_tilde.t() * MX_tilde;
    const arma::mat kH_inv = H.i();
    const arma::vec kbeta_upd = kH_inv * g;
    
    // Trace.
    if (kc_trace > 1) {
      Rcpp::Rcout << "3. Update Fixed Effects" << std::endl;
    }
    
    // Update contribution of the fixed effects.
    b = y_tilde - X_tilde * kbeta_upd - kMy_tilde + MX_tilde * kbeta_upd;
    const arma::vec kD_alpha_upd = b / w_tilde;
    
    // Correction.
    const double kdev_old = dev;
    double rho = 1.0;
    dev = Deviance(beta + rho * kbeta_upd,
                   D_alpha + rho * kD_alpha_upd,
                   ky, kX, kfamily);
    while (dev > kdev_old && rho > kc_rho_tol) {
      // Update rho.
      rho *= 0.5;

      // Compute deviance.
      dev = Deviance(beta + rho * kbeta_upd,
                     D_alpha + rho * kD_alpha_upd,
                     ky, kX, kfamily);
    }
    
    // Update \beta.
    beta += rho * kbeta_upd;
    
    // Trace.
    if (kc_trace > 1) {
      Rcpp::Rcout << "Deviance= " << dev << ", rho= " << rho << std::endl;
    }
    
    // Check termination (Gradient - convergence).
    if (arma::norm(g) < kc_grad_tol) {
      if (kc_trace > 0) {
        Rcpp::Rcout << "Converged - Small Gradient." << std::endl;
      }
      conv = 0;
      b += w_tilde % D_alpha;
      D_alpha += rho * kD_alpha_upd;
      ll = LogLikelihood(beta, D_alpha, ky, kX, kfamily);
      G = MX_tilde.each_col() % kMy_tilde;
      break;
    }
    
    // Check termination (Stepsize - convergence).
    if (arma::norm(kbeta_upd) < kc_step_tol) {
      if (kc_trace > 0) {
        Rcpp::Rcout << "Converged - Small Step Size." << std::endl;
      }
      conv = 2;
      b += w_tilde % D_alpha;
      D_alpha += rho * kD_alpha_upd;
      ll = LogLikelihood(beta, D_alpha, ky, kX, kfamily);
      G = MX_tilde.each_col() % kMy_tilde;
      break;
    }
    
    // Check termination (Deviance - convergence).
    if (std::abs(kdev_old - dev) / dev < kc_dev_tol) {
      if (kc_trace > 0) {
        Rcpp::Rcout << "Converged - Small Relative Change in Deviance." <<
          std::endl;
      }
      conv = 3;
      b += w_tilde % D_alpha;
      D_alpha += rho * kD_alpha_upd;
      ll = LogLikelihood(beta, D_alpha, ky, kX, kfamily);
      G = MX_tilde.each_col() % kMy_tilde;
      break;
    }
    
    // Check termination (Reached maximum number of iterations).
    if (iter == kc_iter_max) {
      if (kc_trace > 0) {
        Rcpp::Rcout << "Reached maximum number of iterations." << std::endl;
      }
      conv = 1;
      b += w_tilde % D_alpha;
      D_alpha += rho * kD_alpha_upd;
      ll = LogLikelihood(beta, D_alpha, ky, kX, kfamily);
      G = MX_tilde.each_col() % kMy_tilde;
      break;
    }
    
    // Update D\alpha.
    D_alpha += rho * kD_alpha_upd;
    
    // Increase 'iter'.
    iter++;
  }
  
  // Return result list.
  return Rcpp::List::create(Rcpp::Named("coefficients") = beta,
                            Rcpp::Named("D.alpha") = D_alpha,
                            Rcpp::Named("maximum") = ll,
                            Rcpp::Named("deviance") = dev,
                            Rcpp::Named("b") = b,
                            Rcpp::Named("w.tilde") = w_tilde,
                            Rcpp::Named("gradient") = g,
                            Rcpp::Named("gradient.cont") = G,
                            Rcpp::Named("Hessian") = H,
                            Rcpp::Named("conv") = conv,
                            Rcpp::Named("iter") = iter);
}