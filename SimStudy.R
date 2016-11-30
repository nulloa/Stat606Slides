library(laGP)
require(MASS)
require(plyr)
require(reshape2)
require(ggplot2)

# Set a seed for repeatable plots
set.seed(12345)


## build up a design with N=~40K locations
n <- 10
x.star <- seq(-n, n, by=1)
X <- as.matrix(expand.grid(x.star, x.star))
#plot(X[,1],X[,2])

calcSigma <- function(X,theta=1){
  Sigma <- matrix(rep(0, nrow(X)*nrow(X)), nrow=nrow(X))
  for (i in 1:nrow(Sigma)) {
    for (j in 1:ncol(Sigma)) {
      Sigma[i,j] <- exp(-(sqrt((X[i,1]-X[j,1])^2+(X[i,2]-X[j,2])^2)^2/theta))
    }
  }
  return(Sigma)
}


## creates the response
getresp <- function(x,theta){
  sigma <- calcSigma(x,theta)
  values <- matrix(rep(0,nrow(x)), ncol=1)
  values <- mvrnorm(1, rep(0, nrow(x)), sigma)
  values <- data.frame(y=values)
}

## local analysis, first pass
Xref <- matrix(c(-2, 2), nrow=TRUE)


## Initialize
iter        <- 100
theta       <- seq(0.01, 10, length=iter)
truepred    <- rep(NA, iter)
nnThetaEst  <- rep(NA, iter)
nnPred      <- rep(NA, iter)
alcThetaEst <- rep(NA, iter)
alcPred     <- rep(NA, iter)

for(i in 1:length(theta)){
  Z0 <- getresp(X, theta[i])$y
  
  truepred[i]    <- Z0[261]
  nngp           <- laGP(Xref, 6, 100, X, Z0, method="nn")
  nnThetaEst[i]  <- nngp$mle$d
  nnPred[i]      <- nngp$mean
  algp           <- laGP(Xref, 6, 100, X, Z0, method="alc", alc.gpu=TRUE)
  alcThetaEst[i] <- algp$mle$d
  alcPred[i]     <- algp$mean
}


SimRes <- cbind(theta,truepred, nnThetaEst, nnPred, alcThetaEst, alcPred)
save(SimRes, file="SimRes.RData")
