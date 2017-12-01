#include <Rcpp.h>
using namespace Rcpp;

// x : sorted vector of values
// w : sorted (w.r.t x ) vector of weights 
// q : sorted vector of quantiles to calculate

// [[Rcpp::export]]
NumericVector weighted_quantile_cpp(NumericVector x, NumericVector w, NumericVector q) {
  double w_sum           = sum(w);
  int x_len              = x.length();
  int q_len              = q.length();
  
  NumericVector w_cumsum    = cumsum(w);
  NumericVector q_threshold = q * w_sum;
  NumericVector q_out(q.length());
  
  int i   = 0;
  int i_q = 0;
  
  for(i = 0; i < x_len; i++) {
    if(w_cumsum[i] > q_threshold[i_q]) {
      q_out[i_q] = x[i];
      i_q++;
      if(i_q == q_len) break;
    }
  }
  
  return q_out;
}
