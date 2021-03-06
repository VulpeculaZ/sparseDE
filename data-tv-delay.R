library(deSolve)
library(CollocInfer)
library(limSolve)
source("./R/tv-delay.R")
source("./R/mDTVSIRfn.R")
source("./R/sparse.R")
source("./R/LS.sparse.R")

## Function to simulate  SIR model:
dtvSIR.gen <- function(t, y, parms){
    if(t < 0){
        lagI1 <- 400
        lagI2 <- 400
    }
    else{
        lagI1 <- lagvalue(t - parms["tau1"], 2)
        lagI2 <- lagvalue(t - parms["tau2"], 2)
    }
    sint <- 8000 * (sin(t / parms["f"] / pi) / 2 + 2)
    dyS <- - (tvtrans(t, kappa)) * (lagI1 + lagI2) * y[1] +  sint
    dyI <- (tvtrans(t, kappa)) * (lagI1 + lagI2) * y[1] - parms["gamma"] * y[2]
    list(c(dyS, dyI))
}

yinit <- c(5000, 800)
tvSIR.pars <- c(10, 1)
tau <- c(4/52 , 5/52)
dtvSIR.pars <- c(tvSIR.pars, tau)
names(dtvSIR.pars) <- c("gamma", "f", "tau1","tau2")
times <- seq(-dtvSIR.pars["tau2"], 5, by = 1/52)
kappa <- c(0.0016, 0.008)
yout <- dede(y = yinit, times = times, func = dtvSIR.gen, parms = dtvSIR.pars, atol = 1e-10)
matplot(yout[,1], yout[,-1], type = "l", lwd = 2, main = "Time Varying SIR Model")

times0 <- knots0 <- times[times >= 0.5]
times.d <- knots.d <- times[times >= 1]
norder = 3
nbasis.d = length(knots.d) + norder - 2
nbasis0 <- length(knots0) + norder - 2
range0  = range(knots0)
range.d <- range(knots.d)
basis0 <- create.bspline.basis(range=range(knots0), nbasis=nbasis0, norder=norder, breaks=knots0)
basis.d <- create.bspline.basis(range=range.d, nbasis=nbasis.d, norder=norder, breaks=knots.d)
fdnames=list(NULL,c('S', 'I'),NULL)
bfdPar0 = fdPar(basis0,lambda=1,int2Lfd(1))
bfdPar.d <- fdPar(basis.d,lambda=1,int2Lfd(1))
## args <- commandArgs(TRUE)
## dataRange <- (1 + 25 * as.numeric(args[1])) : (25 * (as.numeric(args[1]) + 1))
## filename <- paste("nnls-2dadj-sd02-", as.numeric(args[1]),".RData", sep = "")
## initUnif <- runif(100, -1,1)
## initUnifKappa <- runif(100, -0.001,0.001)

set.seed(42)
xout <- c()
xout <- cbind(xout, yout[,2] + rnorm(length(yout[,2]), sd = 200))
xout <- cbind(xout, yout[,3] + rnorm(length(yout[,2]), sd = 200))
xout0 <- xout[times >= 0.5,]
xout.d <- xout[times >= 1,]
DEfd0 <- smooth.basis(knots0, xout0, bfdPar0,fdnames=fdnames)$fd
DEfd.d <- smooth.basis(knots.d, xout.d, bfdPar.d, fdnames=fdnames)$fd
coefs0 <- DEfd0$coefs
colnames(coefs0) = c("S","I")
coefs.d <- DEfd.d$coefs
colnames(coefs.d) = c("S", "I")
lambda = 1000
## Data
dsirData <- matrix(xout.d, ncol =2, dimnames = list(NULL,c("S", "I")))
## Setting initial values

initBeta <- rep(0.1, 10)
initBeta[c(5,6)] <- (initBeta[c(5,6)] + 0.1 )
initBeta <- initBeta / 1.2
initPars <- 10 + runif(1, -1, 1)
names(initPars) <- "gamma"

## initKappa <- rep(0.008, 12) + runif(12, -0.001, 0.001)
## names(initKappa) <- c("k1", "k2", "k3","k4","k5","k6","k7","k8","k9","k10","k11", "k12")

initKappa <- rep(0.012, 4)
names(initKappa) <- c("k1", "k2", "k3","k4")

dde.fit <- Profile.LS.tv.delay(fn = mDTVSIRfn, dsirData, times.d, pars = initPars, beta = initBeta, kappa = initKappa, coefs = coefs.d, basisvals = basis.d, lambda = 1000, in.meth='nlminb', delay = delay, basisvals0 = basis0, coefs0 = coefs0, nbeta = length(initBeta), ndelay = 2, tau = list(seq(0,9/52, by = 1/52)), control.out = list(method = "nnls", maxIter = 15, lambda.sparse = 0))

DEfd.fit <- fd(dde.fit$res$coefs[,2, drop = FALSE], basis.d)
pdf(width = 7, height = 5)
plotfit.fd(dsirData[,2],times.d,DEfd.fit)

DEfd.fit <- fd(dde.fit$res$coefs[,1, drop = FALSE] ,basis.d)
plotfit.fd(y = dsirData[,2] , argvals = times.d, fdobj = DEfd.fit)
dev.off()

## Simplest case

dtvSIR.gen <- function(t, y, parms, kappa){
    if(t < 0){
        lagI1 <- 2000
        ## lagI2 <- 400
    }
    else{
        lagI1 <- lagvalue(t - parms["tau1"], 2)
    }
    sint <- 8000 * (sin(t / parms["f"] / pi) / 2 + 2)
    dyS <- - (tvtrans(t, kappa)) * (lagI1) * y[1] +  sint
    dyI <- (tvtrans(t, kappa)) * (lagI1) * y[1] - parms["gamma"] * y[2]
    list(c(dyS, dyI))
}

yinit <- c(4000, 2000)
tvSIR.pars <- c(10, 1)
tau <- c(8/52)
dtvSIR.pars <- c(tvSIR.pars, tau)
names(dtvSIR.pars) <- c("gamma", "f", "tau1")
times <- seq(-dtvSIR.pars["tau1"], 5, by = 1/52)
kappa <- c(rep(0.005, 3), rep( 0.0025, 3))
yout <- dede(y = yinit, times = times, func = dtvSIR.gen, parms = dtvSIR.pars, atol = 1e-10, kappa = kappa)
matplot(yout[,1], yout[,-1], type = "l", lwd = 2, main = "Time Varying SIR Model")


times0 <- knots0 <- times[times >= 0.5]
times.d <- knots.d <- times[times >= 1]
norder = 3
nbasis.d = length(knots.d) + norder - 2
nbasis0 <- length(knots0) + norder - 2
range0  = range(knots0)
range.d <- range(knots.d)
basis0 <- create.bspline.basis(range=range(knots0), nbasis=nbasis0, norder=norder, breaks=knots0)
basis.d <- create.bspline.basis(range=range.d, nbasis=nbasis.d, norder=norder, breaks=knots.d)
fdnames=list(NULL,c('S', 'I'),NULL)
bfdPar0 = fdPar(basis0,lambda=1,int2Lfd(1))
bfdPar.d <- fdPar(basis.d,lambda=1,int2Lfd(1))
## args <- commandArgs(TRUE)
## dataRange <- (1 + 25 * as.numeric(args[1])) : (25 * (as.numeric(args[1]) + 1))
## filename <- paste("nnls-2dadj-sd02-", as.numeric(args[1]),".RData", sep = "")
## initUnif <- runif(100, -1,1)
## initUnifKappa <- runif(100, -0.001,0.001)


set.seed(42)
xout <- c()
xout <- cbind(xout, yout[,2] + rnorm(length(yout[,2]), sd = 100))
xout <- cbind(xout, yout[,3] + rnorm(length(yout[,2]), sd = 100))
xout0 <- xout[times >= 0.5,]
xout.d <- xout[times >= 1,]
DEfd0 <- smooth.basis(knots0, xout0, bfdPar0,fdnames=fdnames)$fd
DEfd.d <- smooth.basis(knots.d, xout.d, bfdPar.d, fdnames=fdnames)$fd
coefs0 <- DEfd0$coefs
colnames(coefs0) = c("S","I")
coefs.d <- DEfd.d$coefs
colnames(coefs.d) = c("S", "I")
lambda = 1000
## Data
dsirData <- matrix(xout.d, ncol =2, dimnames = list(NULL,c("S", "I")))
## Setting initial values

initBeta <- rep(1/6, 6)
initPars <- 10 + runif(1, -1, 1)
names(initPars) <- "gamma"
initKappa <- kappa
names(initKappa) <- paste("k", 1:length(initKappa), sep = "")

dde.fit <- Profile.LS.tv.delay(fn = mDTVSIRfn, dsirData, times.d, pars = initPars, beta = initBeta, kappa = initKappa, coefs = coefs.d, basisvals = basis.d, lambda = 1000, in.meth='nlminb', delay = delay, basisvals0 = basis0, coefs0 = coefs0, nbeta = length(initBeta), ndelay = 2, tau = list(seq(0, 10/52, by = 2 / 52)), control.out = list(method = "nnls", maxIter = 15, lambda.sparse = 0))

