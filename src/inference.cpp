#include <RcppArmadillo.h>
#include "synthetic.h"
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
// [[Rcpp::export]]
arma::mat inference(arma::mat &y0, arma::mat &y1, arma::mat &x0_scaled,
                    arma::mat &x1_scaled, arma::mat &z0, arma::mat &z1)
{
  int N_ct = x0_scaled.n_cols;
  int T = y0.n_rows;
  arma::mat result(T, N_ct);
  //iterate over units - leave one at a time
  for (int i_unit = 0; i_unit < N_ct; ++i_unit) {
    
    //first generate copies of the input data with the unit left out
    arma::mat x0_scaled_tmp(x0_scaled);
    arma::mat z0_tmp(z0);
    arma::mat y0_tmp(y0);
    
    x0_scaled_tmp.shed_col(i_unit);
    z0_tmp.shed_col(i_unit);
    y0_tmp.shed_col(i_unit);
    
    
    Rcpp::List out = synthetic(y0_tmp, y1, x0_scaled_tmp, x1_scaled, z0_tmp, z1);
    arma::colvec ATT = out["ATT"];
    result.col(i_unit) = ATT;
  }
  
  return result;
  
} 