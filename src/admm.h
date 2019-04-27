#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;

//' quadratic programming with linear equality and box inequality constraints using ADMM
//'@param Dmat : a size (n x n) matrix representing D matrix involved in the quadratic term x. Must be positive semi-definite.
//'@param dvec : a size n vector representing d involved in the linear term of x.
//'@param Amat : a size (m x n) matrix representing A matrix involved in the linear equality constraint of x. 
//'@param bvec : a size m vector representing b involved in the linear equality constraint.
//'@param lb   : a size n vector representing the lower bound of box-inequality constraint of x.
//'@param ub   : a size n vector representing the upper bound of box-inequality constraint of x.
//'@return : a vector containing x, where x is a size N_co vector containing optimal weights for control units.
//'@export
arma::colvec admm (
    arma::mat Dmat, 
    arma::mat dvec,
    arma::mat Amat, 
    arma::mat bvec, 
    arma::mat lb, 
    arma::mat ub
);