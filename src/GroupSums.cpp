#define ARMA_NO_DEBUG
#include <RcppArmadillo.h>


// [[Rcpp::export(name = "groupSums")]]
arma::vec GroupSums(const arma::mat &kM,
                    const arma::vec &kw,
                    const arma::ivec &ka,
                    const arma::ivec &kb) {
  // Auxiliary variables
  const int kn = kM.n_rows;
  const int kp = kM.n_cols;
  const int kl = ka.max() + 1;
  
  // Sort by grouping variable
  arma::ivec a(kn);
  arma::mat M(kn, kp);
  arma::vec w(kn);
  for (int i = 0 ; i < kn ; ++i) {
    const int kb_i = kb(i);
    a(i) = ka(kb_i);
    M.row(i) = kM.row(kb_i);
    w(i) = kw(kb_i);
  }

  // Compute sum of weighted group sums using sorted data
  arma::vec b(kp, arma::fill::zeros);
  int i = 0;
  for (int j = 0 ; j < kl ; ++j) {
    // Compute numerator and denominator
    arma::rowvec num(kp, arma::fill::zeros);
    double denom = 0.0;
    while (i < kn && a(i) == j) {
      num += M.row(i);
      denom += w(i);
      ++i;
    }
    
    // Add weighted group sum
    b += (num / denom).t();
  }
  
  // Return vector
  return b;
}


// [[Rcpp::export(name = "groupSumsSpectral")]]
arma::vec GroupSumsSpectral(const arma::mat &kM,
                            const arma::vec &kv,
                            const arma::vec &kw,
                            const int kL,
                            const arma::ivec &ka,
                            const arma::ivec &kb) {
  // Auxiliary variables
  const int kn = kM.n_rows;
  const int kp = kM.n_cols;
  const int kl = ka.max() + 1;
  
  // Sort by grouping variable
  arma::ivec a(kn);
  arma::mat M(kn, kp);
  arma::vec v(kn);
  arma::vec w(kn);
  for (int i = 0 ; i < kn ; ++i) {
    const int kb_i = kb(i);
    a(i) = ka(kb_i);
    M.row(i) = kM.row(kb_i);
    v(i) = kv(kb_i);
    w(i) = kw(kb_i);
  }
  
  // Compute sum of weighted group sums using sorted data
  arma::vec b(kp, arma::fill::zeros);
  int i = 0;
  for (int j = 0 ; j < kl ; ++j) {
    // Subset individual 'j'
    const int ki1 = i;
    while (i < kn && a(i) == j) ++i;
    const int kTi = i - ki1;
    
    // Compute numerator given a bandwidth 'kL'
    arma::rowvec num(kp, arma::fill::zeros);
    for (int l = 1 ; l <= kL ; ++l) {
      for (int t = l ; t < kTi ; ++t) {
        num += M.row(ki1 + t) * v(ki1 + t - l) * kTi / (kTi - l);
      }
    }
    
    // Compute denominator
    double denom = 0.0;
    for (int t = 0 ; t < kTi ; ++t) {
      denom += w(ki1 + t);
    }
    
    // Add weighted group sum
    b += (num / denom).t();
  }
  
  // Return vector
  return b;
}


// [[Rcpp::export(name = "groupSumsVar")]]
arma::mat GroupSumsVar(const arma::mat &kM,
                       const arma::ivec &ka,
                       const arma::ivec &kb) {
  // Auxiliary variables
  const int kn = kM.n_rows;
  const int kp = kM.n_cols;
  const int kl = ka.max() + 1;
  
  // Sort by grouping variable
  arma::ivec a(kn);
  arma::mat M(kn, kp);
  for (int i = 0 ; i < kn ; ++i) {
    const int kb_i = kb(i);
    a(i) = ka(kb_i);
    M.row(i) = kM.row(kb_i);
  }
  
  // Compute group sum using sorted data
  arma::mat V(kp, kp, arma::fill::zeros);
  int i = 0;
  for (int j = 0 ; j < kl ; ++j) {
    arma::rowvec v(kp, arma::fill::zeros);
    while (i < kn && a(i) == j) {
      v += M.row(i);
      ++i;
    }
    V += v.t() * v;
  }
  
  // Return matrix 
  return V;
}


// [[Rcpp::export(name = "groupSumsCov")]]
arma::mat GroupSumsCov(const arma::mat &kM,
                       const arma::mat &kN,
                       const arma::ivec &ka,
                       const arma::ivec &kb) {
  // Auxiliary variables
  const int kn = kM.n_rows;
  const int kp = kM.n_cols;
  const int kl = ka.max() + 1;
  
  // Sort by grouping variable
  arma::ivec a(kn);
  arma::mat M(kn, kp);
  arma::mat N(kn, kp);
  for (int i = 0 ; i < kn ; ++i) {
    const int kb_i = kb(i);
    a(i) = ka(kb_i);
    M.row(i) = kM.row(kb_i);
    N.row(i) = kN.row(kb_i);
  }
  
  // Compute group sums using sorted data
  arma::mat V(kp, kp, arma::fill::zeros);
  int i = 0;
  for (int j = 0 ; j < kl ; ++j) {
    const int ki1 = i;
    while (i < kn && a(i) == j) ++i;
    const int kTi = i - ki1;
    for (int t = 0 ; t < kTi ; ++t) {
      for (int s = t + 1 ; s < kTi ; ++s) {
        V += M.row(ki1 + t).t() * N.row(ki1 + s);
      }
    }
  }
  
  // Return matrix 
  return V;
}
