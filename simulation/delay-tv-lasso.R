source("./R/timevar.R")
source("./R/tv-delay.R")
source("./R/mDTVSIRfn.R")
source("./R/sparse.R")
source("./R/LS.sparse.R")

library(penalized)
library(CollocInfer)
library(deSolve)
library(nnls)

args <- commandArgs(TRUE)
dataRange <- (1 + 25 * as.numeric(args[1])) : (25 * (as.numeric(args[1]) + 1))
nnls.filename <- paste("nnls-6d-6tv-sd100-", as.numeric(args[1]),".RData", sep = "")
res.filename <- paste("tv-lasso-6d-6tv-sd100-", as.numeric(args[1]),".RData", sep = "")
dataRange <- (1 + 25 * as.numeric(args[1])) : (25 * (as.numeric(args[1]) + 1))

load(nnls.filename)
load("data-tv-1d-sd100.RData")

tau <- 8/52
times <- seq(tau, 5, by = 1/52)
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

xout0 <- data.res[[1]]$xout[times >= 0.5,]
DEfd0 <- smooth.basis(knots0, xout0, bfdPar0,fdnames=fdnames)$fd
coefs0 <- DEfd0$coefs

set.seed(42)
sim.res <- list()
for(i in 1:length(nnls.res)){
    print(dataRange[i])
    initPars <- nnls.res[[i]]$pars
    initBeta <- nnls.res[[i]]$conv$pars.kappa.beta[dim(nnls.res[[i]]$conv$pars.kappa.beta)[1], 8:13]
    initKappa <- nnls.res[[i]]$kappa
    xout0 <- data.res[[dataRange[i]]]$xout[times >= 0.5,]
    xout.d <- data.res[[dataRange[i]]]$xout[times >= 1,]
    DEfd0 <- smooth.basis(knots0, xout0, bfdPar0,fdnames=fdnames)$fd
    coefs0 <- DEfd0$coefs
    nnls.resi <- nnls.res[[i]]
    res.tv.delay <- sparse.tv.delay(fn = mDTVSIRfn, xout.d, times = times.d, pars = initPars, beta = initBeta, kappa = initKappa, basisvals = basis.d, lambda = 1000, in.meth='nlminb',  delay = delay, basisvals0 = basis0, coefs0 = coefs0, control.out = list(method = "fused", maxIter = 10, lambda.sparse = -1), nbeta = length(initBeta), ndelay = 2, tau = list(seq(0, 10/52, by = 2 / 52)), nnls.res = nnls.res[[i]])
    sim.res[[i]] <- res.tv.delay
    save(sim.res, file = res.filename)
}
