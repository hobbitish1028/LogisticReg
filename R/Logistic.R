#'Binomial Logistic regression
#'
#'Fit a generalized linear model via penalized maximum likelihood.
#'The regularization path is computed for the lasso or elasticnet penalty
#'at a grid of values for the regularization parameter lambda.
#'Can deal with all shapes of data,including very large sparse data matrices.
#'
#' @param X input matrix, of dimension n (sample number) by p (variable number);
#' each row is an observation vector.
#' Can be in sparse matrix format
#' (inherit from class "sparseMatrix" as in package Matrix;
#' @param y response variable.
#' Since we aim at binomial distribution, it should be either a factor with two levels,
#'  or a two-column matrix (every row is the probability of two class.)
#' @param max_ite the Maximum number of iterations when we use optimization to estimate the parameeter;
#'  default is 10^5.
#'  @param alpha The elasticnet mixing parameter, with $0\le \alpha \le 1$. The penalty is defined as
#'  $(1-\alpha)/2||\beta||_2^2+\alpha||\beta||_1$.
#'  $\alpha=1$ is the lasso penalty, and $\alpha=0$ the ridge penalty.
#'
#'  @return a data frame which contains the result of this logistic regression.
#'
#'   @param beta the regression result, a length p+1 vector (include the intercept)
#'
#'   @param loss the record of loss in iterations, which can help users to check the convergence
#'
#'   @param Train_Acc the accuracy of the train set, defined as (number of right prediction divided by sample size)


Logistic<-function(X, y, max_ite = 5000, alpha = 1){
  
  m_nosadam<-0
  v_nosadam<-1
  x<-rep(0,p+1)
  X<-cbind(rep(1,dim(X)[1]),X)
  #test_x<-cbind(rep(1,dim(test_x)[1]),test_x)
  alpha_nosadam<-5e-2
  beta_1<-0.9
  #times<-max_ite
  times<-5000
  gamma<-1e-2
  B<-cumsum((1:(1+times))^(-gamma))
  loss<-rep(0,times)
  error<-rep(0,times)
  lambda<-exp(-8)
  
  for(i in 1:times){
    beta_2<-B[i]/B[i+1]
    gradient<-Gradient(X,x,y)
    m_nosadam <- m_nosadam * beta_1 + (1-beta_1) * gradient
    v_nosadam <- v_nosadam * beta_2 + (1-beta_2) *gradient^2
    
    tmp<- abs(gradient)>lambda
    x<- tmp * ( x - alpha_nosadam / sqrt(i) * m_nosadam / sqrt(v_nosadam)      )
    
    loss[i]<-LOSS(X,x,y)
    loss[i]<-LOSS(X,x,y)
    if(i>=10 ){
      error[i]<-var(loss[(i-9):i])
      if(error[i]<1e-3){
        break
      }
    }
  }
  
  predict0<- exp(X%*%x)>=1
  train_acc_nosadam<-mean(predict0==y)
  
  result<-list()
  result$beta<-x
  result$loss<-loss[1:i]
  result$Train_Acc<-train_acc_nosadam
  
  return(result)
  
}


Gradient<-function(X,x,y){
  P0<-exp(X%*%x)
  P<<-P0/(1+P0)
  return( -t(X)%*%(y-P) )
}

LOSS<-function(X,x,y){
  P1<-log(P)
  P0<-log(1-P)
  sum( -P1*(y==1) - P0*(y==0)    )
}

