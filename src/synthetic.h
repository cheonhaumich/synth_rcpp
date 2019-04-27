#include <RcppArmadillo.h>
#include "admm.h"
#include "fn_v.h"
//#include <RInside.h>  
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;

//'This is the wrapping function of 1) estimation, 2) inference with parallel computing, and 3) graph
//'@param y0 : a (T x N_co) matrix of outcomes from control units
//'@param y1 : a (T-dim) vector of outcomes from treated units
//'@param x0_scaled : a (K x N_co) matrix of predictors for control units during pre-treatment times
//'@param x1_scaled : a K-dim vector of predictors for treated units during pre-treatment times
//'@param z0 : a (T0 x N_co) matrix of outcomes for control units during pre-treatment times
//'@param z1 : a (T0-dim) vector of outcomes for treated units during pre-treatment times
//'@return A list containing the following values
//' * ATT : (numeric) a T-dim vector of average treatment effect on treated across time (from the beginning to the end of study time)
//' * T_trajectory : (numeric) a T-dim vector of treatment group outcome trajectory across time
//' * SC_trajectory : (numeric) a T-dim vector of synthetic control outcome trajectory across time
//' * solution.v : (numeric) a (k x k) matrix of weight 
//' * solution.w : (numeric) a Nc-dim vector of synthetic weights for control units
//' @export
Rcpp::List synthetic(
                      arma::mat &y0, 
                      arma::mat &y1,
                      arma::mat &x0_scaled, 
                      arma::mat &x1_scaled, 
                      arma::mat &z0, 
                      arma::mat &z1
);