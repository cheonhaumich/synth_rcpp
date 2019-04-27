#include <RcppArmadillo.h>
#include "admm.h"
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
// [[Rcpp::export]]
arma::colvec admm (
                arma::mat Dmat, 
                arma::mat dvec,
                arma::mat Amat, 
                arma::mat bvec, 
                arma::mat lb, 
                arma::mat ub
                ){
  
  double rho = 1; 
  int max_iter = 10000;

  int n = Dmat.n_cols;
  int m = Amat.n_rows;

  arma::colvec z(n), u(n), x(n);
  z.zeros();
  u.zeros();
  x.zeros();

  arma::mat zero(m,m);
  arma::mat rho_I(n,n);
  zero.zeros();
  rho_I.zeros();
  for ( int i = 0; i < n; i ++ ){
    rho_I(i,i) = rho; 
  } 
  
  
  arma::mat b1, b2, x_mat;
  b1 = join_horiz(Dmat + rho_I, Amat.t());
  b2 = join_horiz(Amat, zero);
  

  x_mat = join_vert(b1, b2);

  arma::mat constant(m + n, 1);
  constant.submat(n, 0, n + m - 1, 0) = bvec;
  
  
  double tol = 1e-8;
  double err = 1;
  for ( int i = 0; i < max_iter; i ++ ){

    // update x and lagrangian nu
    constant.submat(0, 0, n - 1, 0) = rho*(z - u) - dvec;
  
    arma::mat x_tmp = arma::solve(x_mat, constant);

      for ( int j = 0; j < n ; j ++ ) {
        x(j) = x_tmp(j);
      }
  
    // update z
    z = x + u;
    arma::uvec indices;
    indices = arma::find( z < lb );
    z.elem(indices) = lb(indices);
    indices = arma::find( z > ub );
    z.elem(indices) = ub(indices);
    // update u
    u = u+(x-z);
    
    err = arma::norm(x - z, 2)/(arma::norm(x, 2) + arma::norm(z, 2));
    if (err < tol) {break;}
  }
  

  return x;
}
