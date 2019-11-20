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


# devtools::install_github("hobbitish1028/Logistic_Reg")

# Rcpp::sourceCpp('src/LogRegCpp.cpp')
# Logistic<-function(X, y, max_ite = 5000, alpha = 1){
#   
#   m_nosadam<-0
#   v_nosadam<-1
#   p<-dim(X)[2]
#   x<-rep(0,p+1)
#   X<-cbind(rep(1,dim(X)[1]),X)
#   #test_x<-cbind(rep(1,dim(test_x)[1]),test_x)
#   alpha_nosadam<-5e-2
#   beta_1<-0.9
#   #times<-max_ite
#   times<-5000
#   gamma<-1e-2
#   B<-cumsum((1:(1+times))^(-gamma))
#   loss<-rep(0,times)
#   error<-rep(0,times)
#   lambda<-exp(-8)
#   
#   for(i in 1:times){
#     beta_2<-B[i]/B[i+1]
#     gradient<-Gradient(X,x,y)
#     m_nosadam <- m_nosadam * beta_1 + (1-beta_1) * gradient
#     v_nosadam <- v_nosadam * beta_2 + (1-beta_2) *gradient^2
#     
#     tmp<- abs(gradient)>lambda
#     x<- tmp * ( x - alpha_nosadam / sqrt(i) * m_nosadam / sqrt(v_nosadam)      )
#     
#     loss[i]<-LOSS(X,x,y)
#     loss[i]<-LOSS(X,x,y)
#     if(i>=10 ){
#       error[i]<-var(loss[(i-9):i])
#       if(!is.na(error[i]) && error[i]<1e-3){
#         break
#       }
#     }
#   }
#   
#   predict0<- exp(X%*%x)>=1
#   train_acc_nosadam<-mean(predict0==y)
#   
#   result<-list()
#   result$beta<-x
#   result$loss<-loss[1:i]
#   result$Train_Acc<-train_acc_nosadam
#   
#   return(result)
#   
# }
# 
# Gradient<-function(X,x,y){
#   P0<-exp(X%*%x)
#   P<<-P0/(1+P0)
#   return( -t(X)%*%(y-P) )
# }
# 
# Gradient2<-function(X,x,y){
#   P0<-exp(eigenMapMatMult(X, x))
#   P<<-P0/(1+P0)
#   return( -t(X)%*%(y-P) )
# }
# 
# LOSS<-function(X,x,y){
#   P1<-log(P)
#   P0<-log(1-P)
#   sum( -P1*(y==1) - P0*(y==0)    )
# }

Logreg<-function(X,y,maxit = 5000){
  library(Rcpp)
  Rcpp::sourceCpp('src/LogRegCpp.cpp')
  n<-dim(X)[1]
  X<-cbind(rep(1,n),X)
  p<-dim(X)[2]
  
  output<-unique(y)
  yy<- as.numeric(y==output[1])

  ### Use rcpp
  result <- LogRegcpp(X,rep(0,p),yy,maxit = maxit)
  result$loss <- result$loss[result$loss !=0 ]
  result$prediction <- result$P > 0.5
  pred<-rep(output[2],n)
  pred[which(result$P > 0.5)]<- output[1]
  result$prediction <- pred
  result$accuracy <- mean(result$prediction == y)
  result$label <-output
  return(result)
}

My_predict<-function(fit,newx){
  n<-dim(newx)[1]
  X<-cbind(rep(1,n),newx)
  p<-dim(X)[2]
  result <- X%*%fit$x
  tmp <-rep(fit$label[2],n)
  tmp[which(result>0)]<-fit$label[1]
  return( tmp)
}



sigma<-4
set.seed(1)
n<-1e4
p<-1e2
mu1<-rnorm(p)
mu2<-rnorm(p)
X1<-matrix(mu1+rnorm(n*p,0,sigma),n,p,byrow = TRUE)
X2<-matrix(mu2+rnorm(n*p,0,sigma),n,p,byrow = TRUE)
### Train data
X<-rbind(X1,X2)
y<-rep(c(1,0),each=n)
### Test data
test_x<-rbind( matrix(mu1+rnorm(10*n*p,0,sigma),10*n,p,byrow = TRUE), 
               matrix(mu2+rnorm(10*n*p,0,sigma),10*n,p,byrow = TRUE)  )
test_y<-rep(c(1,0),each=10*n)

##Accuracy of the training data and testing data

t1<-proc.time()
fit<-Logreg(X,y,maxit = 200)
proc.time()-t1

##########
# fit<-Logistic(X,y)

plot(51:length(fit$loss),fit$loss[-(1:50)],xlab="iteration" ,ylab = "loss" ,main = "Convergence of the result")



plot(fit$P[c(1:100,n+(1:100))],main="Probability (Prediction)")
my_prediction<-My_predict(fit,newx = test_x)
fit$accuracy
mean(my_prediction == test_y)


## GLMNET Package
library(glmnet)
library(glm2)


### result of glm
### train accuracy
t1<-proc.time()
dat<-as.data.frame(cbind(y,X))
fit0<-glm(y~.,data=dat,family=binomial(link=logit))
proc.time()-t1

(mean((fit0$fitted.values[1:n]>0.5) ) + mean((fit0$fitted.values[n+(1:n)]<0.5) ) )/2

pred_glm <-cbind(rep(1,n),test_x) %*% fit0$coefficients
#plot(pred_glm)
mean( pred_glm*rep(c(1,-1),each=10*n)>0  )

round(fit0$coefficients[1:20]*100)-round(fit$x[1:20]*100)




### result of glmnet 
### train accuracy
t1<-proc.time()
fit0<-glmnet(X,y,family = "binomial")
result2<- predict(fit0, newx = X,  type = "class")
proc.time()-t1
mean(result2==y)
### test accuracy
result2<- predict(fit0, newx = test_x,  type = "class")
mean(result2==test_y)

### result of cv.glmnet
### train accuracy
t1<-proc.time()
fit0<-cv.glmnet(X,y,family = "binomial")
result2<- predict(fit0, newx = X,  type = "class")
proc.time()-t1
mean(result2==y)
### test accuracy
result2<- predict(fit0, newx = test_x,  type = "class")
mean(result2==test_y)

library(dplyr)
library(ggplot2)

class <- rep(c("train","test"),3)
accuracy<-c(0.95585,0.9564,0.9208,0.92072,0.9207,0.95575)
method<- rep(c("My_LogReg","glmnet","cv.glmnet"),each=2)
dat<-data.frame(class,accuracy,method)

ggplot(
  dat %>%
    filter(
      class %in% c("train","test")
    ) %>% group_by(method,class) %>% summarize(
      acc = accuracy
    )
)+ aes(
  x=method,
  y=acc)+
  labs(
    x="Method",
    y="Accuracy"
  )+ geom_col(
    aes(fill=factor(class)),
    position='dodge'
  )+coord_cartesian( ylim = c(0.90, 0.97))+ggtitle("Comparison between three models") +
  theme(plot.title = element_text(hjust = 0.5))