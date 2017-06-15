// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// gpdfit_C
Rcpp::List gpdfit_C(NumericVector x);
RcppExport SEXP loo_gpdfit_C(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(gpdfit_C(x));
    return rcpp_result_gen;
END_RCPP
}
// waic_c
arma::mat waic_c(arma::mat log_like);
RcppExport SEXP loo_waic_c(SEXP log_likeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type log_like(log_likeSEXP);
    rcpp_result_gen = Rcpp::wrap(waic_c(log_like));
    return rcpp_result_gen;
END_RCPP
}
