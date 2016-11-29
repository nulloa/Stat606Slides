library(laGP)
library(ggplot2)

## build up a design with N=~40K locations
x <- seq(-2, 2, by=0.02)
X <- as.matrix(expand.grid(x, x))
Z <- f2d(X)
## local analysis, first pass
Xref <- matrix(c(-1.725, 1.725), nrow=TRUE)

theta <- seq(-2, 2, by=0.2)
for(i in 1:length(theta)){
  nnThetaEst[i]  <- laGP(Xref, 6, 500, X, Z, d=theta[i], method="nn")$mle$d
  nnPred[i]      <- laGP(Xref, 6, 500, X, Z, d=theta[i], method="nn")$mean
  alcThetaEst[i] <- laGP(Xref, 6, 500, X, Z, d=theta[i], method="alc", alc.gpu=TRUE)$mle$d
  alcPred[i]     <- laGP(Xref, 6, 500, X, Z, d=theta[i], method="alc", alc.gpu=TRUE)$mean
}


SimRes <- cbind(nnThetaEst, nnPred, alcThetaEst, alcPred)
save(SimRes, file="SimRes.RData")
