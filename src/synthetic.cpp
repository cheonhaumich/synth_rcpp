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
// [[Rcpp::export]]
Rcpp::List synthetic(
    arma::mat &y0, 
    arma::mat &y1,
    arma::mat &x0_scaled, 
    arma::mat &x1_scaled, 
    arma::mat &z0, 
    arma::mat &z1
){   
  
  Rcpp::Rcout << "start...\n";
  arma::uvec indices;
  
  
  // optimization 1 (V weight)
  // starting value using equal weights
  
  
  
  Rcpp::Rcout << "searching for synthetic control unit... \n";
  
  
  Rcpp::Environment optimxfunction("package:optimx");
  Rcpp::Function optimx = optimxfunction["optimx"];
  //Rcpp::Function collect = optimxfunction["collect.optimx"];
  
  int nrow = x0_scaled.n_rows;
  Rcpp::NumericVector sv(nrow);
  for (int i = 0; i < nrow; i++)
  {
    sv[i] = 1.0/nrow;
  }
  
  Rcpp::List rvgoptim = optimx( Rcpp::_["par"]= sv, 
                                Rcpp::_["fn"] = Rcpp::InternalFunction(&fn_v), 
                                Rcpp::_["method"] = "BFGS", 
                                Rcpp::_["itnmax"] = 10000,
                                Rcpp::_["x0_scaled"] = x0_scaled,
                                Rcpp::_["x1_scaled"] = x1_scaled,
                                Rcpp::_["z0"] = z0,
                                Rcpp::_["z1"] = z1);
  
  arma::colvec sol_v(nrow);
  for (int i = 0; i < nrow; ++i){
    sol_v[i]  = rvgoptim[i];
  }
  
  arma::colvec solution_v;
  
  solution_v = abs(sol_v)/ sum( abs(sol_v) );
  
  Rcpp::Rcout << "V optimization ended \n";
  Rcpp::Rcout << "W optimization begins \n";
  
  
  
  // optimization 2 (W weight)
  // recover w
  arma::mat V(nrow, nrow);
  V.zeros();
  for (int i = 0; i < nrow ; i ++ ){ 
    V(i,i) = solution_v(i);
  }
  
  
  arma::mat H = x0_scaled.t() * V * x0_scaled;
  arma::mat c = -1 * (x1_scaled.t() * V * x0_scaled);
  arma::mat A(1, c.n_cols);
  A.ones();
  arma::vec l(c.n_cols);
  l.zeros();
  arma::vec u(c.n_cols);
  u.ones();
  
  
  
  // quadratic programming
  arma::colvec solution_w;
  arma::mat b_vec_tmp(1,1);
  b_vec_tmp.ones();
  solution_w = admm(H, c.t(), A, b_vec_tmp, l, u); 
  
  Rcpp::Rcout << "W optimization ended \n";
  Rcpp::Rcout << "ATT trajectory is under calculation...\n";
  
  
  // ATT 
  arma::colvec T_trajectory, SC_trajectory, ATT; 
  T_trajectory = y1; 
  SC_trajectory = y0 * solution_w; 
  ATT = T_trajectory - SC_trajectory;
  
  Rcpp::List out = Rcpp::List::create(
    Named("V") = solution_v, 
    Named("W") = solution_w,
    Named("T_trajectory") = T_trajectory,
    Named("SC_trajectory") = SC_trajectory,
    Named("ATT") = ATT
  );
  return out;
}