#include "fn_v.h"
#include "admm.h"
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
// [[Rcpp::export]]
double fn_v(Rcpp::NumericVector &v_vec, 
                 arma::mat x0_scaled,  
                 arma::mat x1_scaled, 
                 arma::mat z0, 
                 arma::mat z1){
  
 arma::vec variables_v(as<arma::vec>(v_vec));
 int number = variables_v.size();
 

 // rescale par
 arma::mat V(number, number); 
 V = V.zeros();
 double sum = 0;
 for ( int i = 0; i < number; i++ ) {
   sum = sum + fabs(variables_v(i)); 
   }

 for ( int i = 0; i < number; i++ ) { V(i,i) = fabs(variables_v(i))/sum; }


 // set up QP
 arma::mat H = x0_scaled.t() * V * x0_scaled; 
 
 arma::mat c = -1 * (x1_scaled.t() * V * x0_scaled);
 int N_ct = c.n_cols;
 arma::mat A(1,N_ct);
 A = A.ones();
 arma::vec l( N_ct );
 l = l.zeros();
 arma::vec u( N_ct );
 u = u.ones();
 arma::vec b(1);  
 b = b.ones();
  
  // quadratic programming
 arma::colvec solution_w; 
 solution_w = admm(H, c.t(), A, b, l, u);


 
 // compute loses
 arma::mat loss_w(1,1);
 loss_w = ( x1_scaled - x0_scaled * solution_w ).t() * V * ( x1_scaled - x0_scaled * solution_w ); 
 arma::mat loss_v(1,1);

 loss_v = ( z1 - z0 * solution_w ).t() * ( z1 - z0 * solution_w );
 loss_v = loss_v / z0.n_elem;

 return loss_v(0);
}