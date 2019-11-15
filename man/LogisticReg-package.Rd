\name{LogisticReg}
\alias{LogisticReg}
\alias{LogisticReg}
\docType{package}
\title{
\packageTitle{LogisticReg}
}
\description{
\packageDescription{LogisticReg}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{LogisticReg}
\packageIndices{LogisticReg}
 Logreg(X,y,maxit = 5000)
 X is the n by p data matrix, and y is binary n by 1 response vector. In the optimization part, the maximal iteration steps is setted as 5000 as default.
}
\author{
\packageAuthor{Chang Wang}

Maintainer: \packageMaintainer{Chang Wang}
}

\keyword{ package, logistic regression }

\examples{

sigma<-5
set.seed(1)
n<-1e4
p<-1e2
mu1<-rnorm(p)
mu2<-rnorm(p)
X1<-matrix(mu1+rnorm(n*p,0,sigma),n,p,byrow = TRUE)
X2<-matrix(mu2+rnorm(n*p,0,sigma),n,p,byrow = TRUE)
X<-rbind(X1,X2)
y<-rep(c("a","b"),each=n)
fit<-Logreg(X,y)
}
