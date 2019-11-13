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

  return(result)
  
}


options(warnings = -1)
library(Rcpp)
library(RcppEigen)

sourceCpp("./src/Matrix_multiply.cpp")
sourceCpp("./src/rcppeigen_hello_world.cpp")

aaa(matrix(1:6,3,2),1:2)

bbb(matrix(1:6,3,2),3:1)

eigenadd(0,matrix(1:4,2,2))

multiply(matrix(1:4,2,2),1:2)

rcppeigen_hello_world()
rcppeigen_bothproducts(1:3)

sigma<-4

set.seed(1)
n<-1e2
p<-1e2
mu1<-rnorm(p)
mu2<-rnorm(p)
X1<-matrix(mu1+rnorm(n*p,0,sigma),n,p,byrow = TRUE)
X2<-matrix(mu2+rnorm(n*p,0,sigma),n,p,byrow = TRUE)
X<-rbind(X1,X2)

y<-rep(c(1,0),each=n)

test_x<-rbind( matrix(mu1+rnorm(n*p,0,sigma),n,p,byrow = TRUE), matrix(mu2+rnorm(n*p,0,sigma),n,p,byrow = TRUE)  )
test_y<-rep(c(1,0),each=n)


t1<-proc.time()
result1<-LogRegcpp(X,rep(0,p),y,maxit=5000)
proc.time()-t1


t1<-proc.time()
result1<-Logistic(X,y)
proc.time()-t1

plot(result1$loss)


