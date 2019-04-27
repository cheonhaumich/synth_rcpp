#include <RcppArmadillo.h>


// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;

//'This is the optimizing function (minimization) for V
//'
//'@param v_vec : a k-dim vector of initial values for V
//'@param x0_scaled : a (k x Nc) matrix of scaled predictors for control group during pre-treatment times
//'@param x1_scaled : a (k x 1) scaled predictors for treated group during pre-treatment times
//'@param z0 : a (T0 x Nc) matrix of outcomes for control group during pre-treatment times
//'@param z1 : a (T0 x 1)  matrix of outcomes for treated group during pre-treatment times
//'@return a scalar of V loss during pre-treatment times : t(X1-X0%*%W) %*% V %*% (X1-X0%*%W)
//'@export
double fn_v(Rcpp::NumericVector &v_vec, 
            arma::mat x0_scaled,  
            arma::mat x1_scaled, 
            arma::mat z0, 
            arma::mat z1);