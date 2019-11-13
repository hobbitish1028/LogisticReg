#include <RcppEigen.h>
#include <math.h>

// [[Rcpp::depends( RcppEigen)]]

// [[Rcpp::export]]
SEXP eigenMatMult(Eigen::MatrixXd A, Eigen::MatrixXd B){
  Eigen::MatrixXd C = A * B;
  
  return Rcpp::wrap(C);
}

// [[Rcpp::export]]
SEXP eigenadd(double A, Eigen::MatrixXd B){
  Eigen::MatrixXd C = B;

  return Rcpp::wrap(C);
}

// [[Rcpp::export]]
SEXP eigenMapMatMult(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B){
  Eigen::MatrixXd C = A * B;
  
  return Rcpp::wrap(C);
}

// [[Rcpp::export]]
Rcpp::NumericVector aaa( Rcpp::NumericMatrix A, Rcpp::NumericVector B){
  return( transpose(A) * B);
}

// [[Rcpp::export]]
Rcpp::NumericVector bbb( Rcpp::NumericVector P0){
  Rcpp::NumericVector PP;
  PP = exp(P0);
  Rcpp::NumericVector P = PP / (1.0 + PP);
  return  P ;
}

// [[Rcpp::export]]
Rcpp:: NumericVector multi(Rcpp:: NumericMatrix X, Rcpp:: NumericVector beta, int nsites, int p){
  
  Rcpp:: NumericVector result(nsites);
  for(int j = 0; j < nsites; j++)
  {
    int temp = 0;
    
    for(int l = 0; l < p; l++) temp = temp + X(j,l) * beta[l];
    
    result[j] = temp;
  }
  
  return result;
}