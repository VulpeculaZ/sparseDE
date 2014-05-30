source("./R/tv-delay.R")
source("./R/mDTVSIRfn.R")
source("./R/sparse.R")
source("./R/LS.sparse.R")
library(CollocInfer)
library(limSolve)
load("data-tv-1d-sd100.RData")

times <- seq(-dtvSIR.pars["tau1"], 5, by = 1/52)
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
args <- commandArgs(TRUE)
dataRange <- (1 + 25 * as.numeric(args[1])) : (25 * (as.numeric(args[1]) + 1))
filename <- paste("nnls-6d-6tv-sd100-", as.numeric(args[1]),".RData", sep = "")

begTime <- Sys.time()
nnls.res <- list()
for(i in 1:length(dataRange)){
    print(dataRange[i])
    ## points(times, xout)
    xout0 <- data.res[[dataRange[i]]]$xout[times >= 0.5,]
    xout.d <- data.res[[dataRange[i]]]$xout[times >= 1,]
    DEfd0 <- smooth.basis(knots0, xout0, bfdPar0,fdnames=fdnames)$fd
    DEfd.d <- smooth.basis(knots.d, xout.d, bfdPar.d, fdnames=fdnames)$fd
    ## temp.fit <- eval.fd(times.d, DEfd.d)
    ## par(ask=FALSE)
    ## plotfit.fd(xout[times >=0,], times0, DEfd0)
    ## plotfit.fd(xout[times >=5,], times.d, DEfd.d)
    ## extract the coefficients and assign variable names
    coefs0 <- DEfd0$coefs
    colnames(coefs0) = c("S","I")
    coefs.d <- DEfd.d$coefs
    colnames(coefs.d) = c("S", "I")
    ##  Set a value for lambda
    lambda = 1000
    ## Data
    dsirData <- matrix(xout.d, ncol =2, dimnames = list(NULL,c("S", "I")))
    ## Setting initial values
    initPars <- data.res[[dataRange[i]]]$initPars
    initBeta <- data.res[[dataRange[i]]]$initBeta
    initKappa <- data.res[[dataRange[i]]]$initKappa
    dde.fit <- Profile.LS.tv.delay(fn = mDTVSIRfn, dsirData, times.d, pars = initPars, beta = initBeta, kappa = initKappa, coefs = coefs.d, basisvals = basis.d, lambda = 1000, in.meth='nlminb', delay = delay, basisvals0 = basis0, coefs0 = coefs0, nbeta = length(initBeta), ndelay = 2, tau = list(seq(0, 10/52, by = 2 / 52)), control.out = list(method = "nnls", maxIter = 15, lambda.sparse = 0))
    nnls.res[[i]] <- dde.fit$res
    save(nnls.res, file =filename)
}

runTime <- Sys.time() - begTime
print(runTime)
save(nnls.res,  runTime, file = filename)
