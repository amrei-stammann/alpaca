#define ARMA_NO_DEBUG
#include <RcppArmadillo.h>


// [[Rcpp::export(name = "groupSums")]]
arma::vec GroupSums(const arma::mat& M,
                    const arma::vec& w,
                    const Rcpp::List& jlist) {
  // Auxiliary variables (fixed)
  const int J = jlist.size();
  const int P = M.n_cols;
  
  // Auxiliary variables (storage)
  double denom;
  int j, p, t, T;
  arma::vec b(P);
  arma::vec num(P);
  
  // Compute sum of weighted group sums
  b.zeros();
  for (j = 0 ; j < J ; ++j) {
    // Subset j-th group
    Rcpp::IntegerVector indexes = jlist[j];
    T = indexes.size();
    
    // Compute numerator of the weighted group sum
    num.zeros();
    for (p = 0 ; p < P ; ++p) {
      for (t = 0 ; t < T ; ++t) {
        num(p) += M(indexes[t], p);
      }
    }
    
    // Compute denominator of the weighted group sum
    denom = 0.0;
    for (t = 0 ; t < T ; ++t) {
      denom += w(indexes[t]);
    }
    
    // Add weighted group sum
    b += num / denom;
  }
  
  // Return vector
  return b;
}


// [[Rcpp::export(name = "groupSumsSpectral")]]
arma::vec GroupSumsSpectral(const arma::mat& M,
                            const arma::vec& v,
                            const arma::vec& w,
                            const int L,
                            const Rcpp::List& jlist) {
  // Auxiliary variables (fixed)
  const int J = jlist.size();
  const int P = M.n_cols;
  
  // Auxiliary variables (storage)
  double denom;
  int j, l, p, t, T;
  arma::vec b(P);
  arma::vec num(P);
  
  // Compute sum of weighted group sums
  b.zeros();
  for (j = 0 ; j < J ; ++j) {
    // Subset j-th group
    Rcpp::IntegerVector indexes = jlist[j];
    T = indexes.size();
    
    // Compute numerator of the weighted group sum given bandwidth 'L'
    num.zeros();
    for (p = 0 ; p < P ; ++p) {
      for (l = 1 ; l <= L ; ++l) {
        for (t = l ; t < T ; ++t) {
          num(p) += M(indexes[t], p) * v(indexes[t - l]) * T / (T - l);
        }
      }
    }
    
    // Compute denominator of the weighted group sum
    denom = 0.0;
    for (t = 0 ; t < T ; ++t) {
      denom += w(indexes[t]);
    }
    
    // Add weighted group sum
    b += num / denom;
  }
  
  // Return vector
  return b;
}


// [[Rcpp::export(name = "groupSumsVar")]]
arma::mat GroupSumsVar(const arma::mat& M,
                       const Rcpp::List& jlist) {
  // Auxiliary variables (fixed)
  const int J = jlist.size();
  const int P = M.n_cols;
  
  // Auxiliary variables (storage)
  int j, p, t, T;
  arma::vec v(P);
  arma::mat V(P, P);
  
  // Compute covariance matrix
  V.zeros();
  for (j = 0 ; j < J ; ++j) {
    // Subset j-th group
    Rcpp::IntegerVector indexes = jlist[j];
    T = indexes.size();
    
    // Compute group sum
    v.zeros();
    for (p = 0 ; p < P ; ++p) {
      for (t = 0 ; t < T ; ++t) {
        v(p) += M(indexes[t], p);
      }
    }
    
    // Add to covariance matrix
    V += v * v.t();
  }
  
  // Return matrix 
  return V;
}


// [[Rcpp::export(name = "groupSumsCov")]]
arma::mat GroupSumsCov(const arma::mat& M,
                       const arma::mat& N,
                       const Rcpp::List& jlist) {
  // Auxiliary variables (fixed)
  const int J = jlist.size();
  const int P = M.n_cols;
  
  // Auxiliary variables (storage)
  int j, p, q, t, s, T;
  arma::mat V(P, P);
  
  // Compute covariance matrix
  V.zeros();
  for (j = 0 ; j < J ; ++j) {
    // Subset j-th group
    Rcpp::IntegerVector indexes = jlist[j];
    T = indexes.size();
    
    // Add to covariance matrix
    for (p = 0 ; p < P ; ++p) {
      for (q = 0 ; q < P ; ++q) {
        for (t = 0 ; t < T ; ++t) {
          for (s = t + 1 ; s < T ; ++s) {
            V(q, p) += M(indexes[t], q) * N(indexes[s], p);
          }
        }
      }
    }
  }
  
  // Return matrix 
  return V;
}