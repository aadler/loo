// [[Rcpp::depends(RcppArmadillo)]]
# include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
arma::mat waic_c(arma::mat log_like){

  arma::mat out(log_like.n_cols, 2);

  for(arma::uword ii = 0; ii<log_like.n_cols; ++ii )
  {
    out(ii,0) = var(log_like.col(ii));
    out(ii,1) = log(mean(exp(log_like.col(ii))));
  }


  return out;
}
