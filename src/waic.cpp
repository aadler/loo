// [[Rcpp::depends(RcppArmadillo)]]
# include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
Rcpp::List waic_c(arma::mat log_like){


  arma::mat::iterator it_end = log_like.end_col(0);
  arma::mat::iterator it_start = log_like.begin_col(log_like.n_cols - 1);

  arma::vec pwaic(log_like.n_cols);
  arma::vec lppd(log_like.n_cols);

  for(arma::uword ii = 0; ii<log_like.n_cols; ++ii )
  {
    pwaic(ii) = var(log_like.col(ii));
    lppd(ii) = log(mean(exp(log_like.col(ii))));
  }


  return Rcpp::List::create((Rcpp::List::create(Rcpp::Named("LPPD") = lppd,
                                                                Rcpp::Named("PWAIC") = pwaic)));
}
