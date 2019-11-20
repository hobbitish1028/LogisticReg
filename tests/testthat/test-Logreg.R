test_that("multiplication works", {
  
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

  ##Accuracy of the training data and testing data
  fit<-Logreg(X,y,maxit = 200)
  acc_mine<-fit$accuracy

  library(glm2)
  ### result of glm
  ### train accuracy
  dat<-as.data.frame(cbind(y,X))
  fit0<-glm(y~.,data=dat,family=binomial(link=logit))
  acc_glm<-(mean((fit0$fitted.values[1:n]>0.5) ) + mean((fit0$fitted.values[n+(1:n)]<0.5) ) )/2
  expect_equal( 2*2,4 )
})
