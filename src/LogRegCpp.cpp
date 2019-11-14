#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;


Rcpp:: NumericVector multi(Rcpp:: NumericMatrix X, Rcpp:: NumericVector beta, int nsites, int p){
  
  Rcpp:: NumericVector result(nsites);
  for(int j = 0; j < nsites; j++)
  {
    double temp = 0;
    
    for(int l = 0; l < p; l++) temp = temp + X(j,l) * beta[l];
    
    result[j] = temp;
  }
  
  return result;
}


// [[Rcpp::export]]
List LogRegcpp(NumericMatrix X, NumericVector x, NumericVector y ,int maxit){
  
  int p = x.size();
  int n = y.size();
  
  NumericVector m(p);
  NumericVector v(p);
  for(int i = 0; i<p ; i++){
    m[i] = 0;
    v[i] = 1;
  } 
  double alpha = 5e-2;
  double beta_1 =0.9;
  double beta_2 = 0.9;
  int times = maxit;
  double gamma = 1e-2;
  
  NumericVector B(times+100); 
  B[0] = 1;
  for(int i=1;i<times+100;i++){
    B[i] = B[i-1] + pow( i,-gamma);  
  }
  
  NumericVector loss(times);
  NumericVector error(times);
  
  
  NumericVector P0(n);
  NumericVector P(n);
  NumericVector gradient(p);
  
  for(int i=1; i < times ; i++){
    beta_2 = B[i+50] / B[i+51];
    NumericMatrix A1 = X;
    NumericVector B1 = x;
    NumericVector P0 =  multi(X ,x,n,p);
    
    for(int i=0;i<n;i++){
      if(P0[i] < -10){
        P0[i] = -10;
      }
      if(P0[i] > 10){
        P0[i] = 10;
      }
    }
    P0 = exp(-P0);
    P = 1.0/ (1.0 + P0);
    
    NumericMatrix A2 = transpose(X);
    NumericVector B2 = P - y;
    //NumericVector gradient =  A2 * B2; 
    NumericVector gradient =  multi(A2 ,B2,p,n);
    
    m = m * beta_1 + (1-beta_1) * gradient;
    v = v * beta_2 + (1-beta_2) * pow(gradient,2.0);
    
    x =   x - alpha / sqrt(i) * m / sqrt(v);
    
    
    NumericVector P1 = log(P);
    NumericVector P2 = log(1-P);
    double loss0 = 0;
    for( int j =0;j<n;j++){
      if(y[j]==0){
        loss0  -= P2[j];
      }else{
        loss0  -= P1[j];
      }
    }
    loss[i] = loss0;  // Loss function is the log likelihood
    
    if( i>20 ){
      double mu = 0; 
      for(int j = i-10; j < i; j++){
        mu += loss[j]; 
      }
      mu = mu /10;
      double s = 0;
      for(int j = i-10; j < i; j++){
        s += pow(loss[i] - mu, 2.0);
      }  
      s = s/10;
      // if the loss converge, the break
      if(s < 0.001){
        break;
      }
      
    }
    
  }
  
  return Rcpp::List::create(Rcpp::Named("x") = x,
                            Rcpp::Named("beta_2") = beta_2,
                            Rcpp::Named("P") = P,
                            Rcpp::Named("m") = m,
                            Rcpp::Named("v") = v,
                            Rcpp::Named("loss") = loss);
}


