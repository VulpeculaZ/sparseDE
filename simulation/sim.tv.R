source("./sources/timevar.R")
## source("./sources/LS.sparse.R")
## source("./sources/DSIRfnSparse.R")
source("./sources/poslasso.R")

library(penalized)
library(CollocInfer)
library(deSolve)
library(nnls)

## Function to simulate DSIR
## Function to simulate  SIR model:

tvSIR.gen <- function(t, y, parms){
    sint <- 2000 * (sin(t / parms["f"]) / 2 + 2)
    dyS <- - (tvtrans(t, kappa)) * y[2] * y[1] +  sint
    dyI <- (tvtrans(t, kappa)) * y[2] * y[1] - parms["gamma"] * y[2]
    list(c(dyS, dyI))
}


yinit <- c(3000, 400)
times <- seq(0, 5, by = 1/52)
tvSIR.pars <- c(20, 1)
kappa <- rep(0.01, 12)
kappa[10:12] <- 0.006
names(tvSIR.pars) <- c("gamma", "f")
yout <- dede(y = yinit, times = times, func = tvSIR.gen, parms = tvSIR.pars, atol = 1e-10)
matplot(yout[,1], yout[,-1], type = "l", lwd = 2, main = "Delayed SIR Model")


knots <- times
norder = 3
nbasis = length(knots) + norder - 2
range  = range(knots)
basis <- create.bspline.basis(range=range(knots), nbasis=nbasis, norder=norder, breaks=knots)
fdnames=list(NULL,c('S', 'I'),NULL)
bfdPar <- fdPar(basis,lambda=1,int2Lfd(1))
initUnif <- runif(100, -1,1)

set.seed(42)
sim.res <- list()
for(i in 1:100){
    xout <- c()
    xout <- cbind(xout, yout[,2] + rnorm(length(yout[,2]), sd = 0.01))
    xout <- cbind(xout, yout[,3] + rnorm(length(yout[,2]), sd = 0.01))
    ## points(times, xout)
    DEfd <- smooth.basis(knots, xout, bfdPar,fdnames=fdnames)$fd
    ## temp.fit <- eval.fd(times.d, DEfd.d)
    ## par(ask=FALSE)
    ## plotfit.fd(xout[times >=0,], times0, DEfd0)
    ## plotfit.fd(xout[times >=5,], times.d, DEfd.d)
    ## extract the coefficients and assign variable names
    coefs <- DEfd$coefs
    colnames(coefs) = c("S","I")
    ##  Set a value for lambda
    lambda = 1
    ## Data
    tvData <- matrix(xout, ncol =2, dimnames = list(NULL,c("S", "I")))
    ## Setting initial values
    initPars <- 20 + initUnif[i]
    names(initPars) <- c("gamma")
    initKappa <- rep(0.01, 12)
    names(initKappa) <- c("k1", "k2", "k3","k4","k5","k6","k7","k8","k9","k10","k11", "k12")
    tv.fit <- Profile.LS.tv(tvDSIRfn, tvData, times, pars = initPars, kappa = initKappa, coefs = coefs, basisvals = basis, lambda = 1, in.meth='nlminb', control.out = list(method = "nnls", maxIter = 10, lambda.sparse = 0))
    sim.res[[i]] <- tv.fit
    save(sim.res, initPars, initKappa,  file ="sim.tv01.RData")
}
