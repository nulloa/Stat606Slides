library(laGP)
require(MASS)
require(plyr)
require(reshape2)
require(ggplot2)

# Set a seed for repeatable plots
set.seed(12345)


## build up a design with N=~40K locations
n <- 10
x.star <- seq(-n, n, length=50)
X <- as.matrix(expand.grid(x.star, x.star))

calcSigma <- function(X1,X2,l=1){
  Sigma <- matrix(rep(0, length(X1)*length(X2)), nrow=length(X1))
  for (i in 1:nrow(Sigma)) {
    for (j in 1:ncol(Sigma)) {
      Sigma[i,j] <- exp(-0.5*(abs(X1[i]-X2[j])/l)^2)
    }
  }
  return(Sigma)
}

## creates the response
getresp <- function(x,theta){
  sigma <- calcSigma(x,x)
  values <- matrix(rep(0,length(x)), ncol=1)
  values <- mvrnorm(1, rep(theta, length(x)), sigma)
  values <- data.frame(x=x, y=values)
}

## local analysis, first pass
Xref <- matrix(c(-1.5, 1.83), nrow=TRUE)


## Initialize
iter        <- 100
theta       <- seq(-1, 1, length=iter)
truepred    <- rep(NA, iter)
nnThetaEst  <- rep(NA, iter)
nnPred      <- rep(NA, iter)
alcThetaEst <- rep(NA, iter)
alcPred     <- rep(NA, iter)

for(i in 1:length(theta)){
  Z1 <- getresp(X[,1], theta[i])$y
  Z2 <- getresp(X[,2], theta[i])$y
  Z0 <- Z1*Z2
  
  truepred[i]    <- Z0[579]
  nngp           <- laGP(Xref, 6, 100, X, Z0, method="nn")
  nnThetaEst[i]  <- nngp$mle$d
  nnPred[i]      <- nngp$mean
  algp           <- laGP(Xref, 6, 100, X, Z0, method="alc", alc.gpu=TRUE)
  alcThetaEst[i] <- algp$mle$d
  alcPred[i]     <- algp$mean
}


SimRes <- cbind(theta,truepred, nnThetaEst, nnPred, alcThetaEst, alcPred)
save(SimRes, file="SimRes.RData")
