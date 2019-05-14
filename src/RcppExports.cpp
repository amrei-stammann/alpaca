// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// CenterVariables
arma::mat CenterVariables(const arma::mat& kV, const arma::vec& kw, const arma::imat& kA, const arma::imat& kB, const arma::ivec& klvls_k, const double ktol);
RcppExport SEXP _alpaca_CenterVariables(SEXP kVSEXP, SEXP kwSEXP, SEXP kASEXP, SEXP kBSEXP, SEXP klvls_kSEXP, SEXP ktolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type kV(kVSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type kw(kwSEXP);
    Rcpp::traits::input_parameter< const arma::imat& >::type kA(kASEXP);
    Rcpp::traits::input_parameter< const arma::imat& >::type kB(kBSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type klvls_k(klvls_kSEXP);
    Rcpp::traits::input_parameter< const double >::type ktol(ktolSEXP);
    rcpp_result_gen = Rcpp::wrap(CenterVariables(kV, kw, kA, kB, klvls_k, ktol));
    return rcpp_result_gen;
END_RCPP
}
// GetAlpha
arma::vec GetAlpha(const arma::vec& kpi, const arma::ivec& klvls_k, const arma::imat& kA, const arma::imat& kB, const double ktol);
RcppExport SEXP _alpaca_GetAlpha(SEXP kpiSEXP, SEXP klvls_kSEXP, SEXP kASEXP, SEXP kBSEXP, SEXP ktolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type kpi(kpiSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type klvls_k(klvls_kSEXP);
    Rcpp::traits::input_parameter< const arma::imat& >::type kA(kASEXP);
    Rcpp::traits::input_parameter< const arma::imat& >::type kB(kBSEXP);
    Rcpp::traits::input_parameter< const double >::type ktol(ktolSEXP);
    rcpp_result_gen = Rcpp::wrap(GetAlpha(kpi, klvls_k, kA, kB, ktol));
    return rcpp_result_gen;
END_RCPP
}
// GroupSums
arma::vec GroupSums(const arma::mat& kM, const arma::vec& kw, const arma::ivec& ka, const arma::ivec& kb);
RcppExport SEXP _alpaca_GroupSums(SEXP kMSEXP, SEXP kwSEXP, SEXP kaSEXP, SEXP kbSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type kM(kMSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type kw(kwSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type ka(kaSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type kb(kbSEXP);
    rcpp_result_gen = Rcpp::wrap(GroupSums(kM, kw, ka, kb));
    return rcpp_result_gen;
END_RCPP
}
// GroupSumsSpectral
arma::vec GroupSumsSpectral(const arma::mat& kM, const arma::vec& kv, const arma::vec& kw, const int kL, const arma::ivec& ka, const arma::ivec& kb);
RcppExport SEXP _alpaca_GroupSumsSpectral(SEXP kMSEXP, SEXP kvSEXP, SEXP kwSEXP, SEXP kLSEXP, SEXP kaSEXP, SEXP kbSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type kM(kMSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type kv(kvSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type kw(kwSEXP);
    Rcpp::traits::input_parameter< const int >::type kL(kLSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type ka(kaSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type kb(kbSEXP);
    rcpp_result_gen = Rcpp::wrap(GroupSumsSpectral(kM, kv, kw, kL, ka, kb));
    return rcpp_result_gen;
END_RCPP
}
// GroupSumsVar
arma::mat GroupSumsVar(const arma::mat& kM, const arma::ivec& ka, const arma::ivec& kb);
RcppExport SEXP _alpaca_GroupSumsVar(SEXP kMSEXP, SEXP kaSEXP, SEXP kbSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type kM(kMSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type ka(kaSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type kb(kbSEXP);
    rcpp_result_gen = Rcpp::wrap(GroupSumsVar(kM, ka, kb));
    return rcpp_result_gen;
END_RCPP
}
// GroupSumsCov
arma::mat GroupSumsCov(const arma::mat& kM, const arma::mat& kN, const arma::ivec& ka, const arma::ivec& kb);
RcppExport SEXP _alpaca_GroupSumsCov(SEXP kMSEXP, SEXP kNSEXP, SEXP kaSEXP, SEXP kbSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type kM(kMSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type kN(kNSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type ka(kaSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type kb(kbSEXP);
    rcpp_result_gen = Rcpp::wrap(GroupSumsCov(kM, kN, ka, kb));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_alpaca_CenterVariables", (DL_FUNC) &_alpaca_CenterVariables, 6},
    {"_alpaca_GetAlpha", (DL_FUNC) &_alpaca_GetAlpha, 5},
    {"_alpaca_GroupSums", (DL_FUNC) &_alpaca_GroupSums, 4},
    {"_alpaca_GroupSumsSpectral", (DL_FUNC) &_alpaca_GroupSumsSpectral, 6},
    {"_alpaca_GroupSumsVar", (DL_FUNC) &_alpaca_GroupSumsVar, 3},
    {"_alpaca_GroupSumsCov", (DL_FUNC) &_alpaca_GroupSumsCov, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_alpaca(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
