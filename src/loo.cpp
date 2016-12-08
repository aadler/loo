#include <Rcpp.h>
using namespace Rcpp;

NumericVector lx_C(NumericVector b, NumericVector y){
  int M = b.size();
  int N = y.size();
  NumericVector a = clone(b);
  NumericVector x = clone(y);
  NumericVector k(M);
  NumericVector res(M);
  for(int i = 0; i < M; ++i)
    a[i] = -a[i];
  double tmp;
  for(int i = 0; i < M; ++i) {
    tmp = 0.0;
    for(int j = 0; j < N; ++j)
      tmp += log1p(a[i] * x[j]);
    k[i] = tmp / N;
  }
  for(int i = 0; i < M; ++i)
    res[i] = log(a[i] / k[i]) - k[i] - 1.0;
  return res;
}

// [[Rcpp::export]]
Rcpp::List gpdfit_C(NumericVector x) {
  int N = x.size();
  NumericVector xx = clone(x);
  std::sort(xx.begin(), xx.end());
  double prior = 3.0;
  int M = 80 + floor(sqrt(N));
  IntegerVector mseq(M);
  for(int i = 0; i < M;)
    mseq[i] = ++i;
  NumericVector sM(M);
  for(int i = 0; i < M; ++i)
    sM[i] = 1.0 - sqrt(M / (mseq[i] - 0.5));
  int Nflr = floor(0.25 * N + 0.5);
  NumericVector b(M);
  double prefix = 1.0 / xx[N - 1];
  double denominator = prior * xx[Nflr - 1];
  for(int i = 0; i < M; ++i)
    b[i] = prefix + sM[i] / denominator;
  NumericVector l(M);
  l = lx_C(b, xx);
  for(int i = 0; i < M; ++i)
    l[i] *= N;
  NumericVector w(M);
  for(int j = 0; j < M; ++j) {
    double tmp = 0.0;
    for(int i = 0; i < M; ++i)
      tmp += exp(l[i] - l[j]);
    w[j] = 1 / tmp;
  }
  double bdotw = 0.0;
  for(int i = 0; i < M; ++i)
    bdotw += b[i] * w[i];
  double tmp = 0.0;
  for(int i = 0; i < N; ++i)
    tmp += log1p(-bdotw * x[i]);
  double k = tmp / N;
  double sigma = -k / bdotw;
  return Rcpp::List::create(Rcpp::Named("k") = k,
                            Rcpp::Named("sigma") = sigma);
}
