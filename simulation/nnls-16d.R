source("./R/DSIRfnSparse.R")

library(spam)
library(gpDDE)
load("data-2dadj-sd02.RData")

times <- seq(-DSIR.pars["tau2"], 25, by = 0.1)
times0 <- knots0 <- times[times >= 0]
times.d <- knots.d <- times[times >= 5]
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
filename <- paste("nnls-2dadj-sd02-true-", as.numeric(args[1]),".RData", sep = "")

begTime <- Sys.time()
set.seed(42)
sim.Covar <- list()
nnls.res <- list()
for(i in 1:length(dataRange)){
    print(dataRange[i])
    ## points(times, xout)
    xout0 <- data.res[[dataRange[i]]]$xout[times >= 0,]
    xout.d <- data.res[[dataRange[i]]]$xout[times >= 5,]
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
    initPars <- 0.5
    names(initPars) <- c("gamma")
    initBeta <- c(2, 2)
    dde.fit <- Profile.LS.DDE(DSIRfn.sparse, dsirData, times.d, pars = initPars, beta = initBeta, coefs = coefs.d, basisvals = basis.d, lambda = 1000, in.meth='nlminb', basisvals0 = basis0, coefs0 = coefs0, nbeta = length(initBeta), ndelay = 2, tau = list(c(2, 7/3)), control.out = list(method = "nnls.old", maxIter = 20, echo = TRUE, tol = 1e-10))
    Covar <- NA
    try(Covar <- ProfileSSE.covariance.DDE(fn = DSIRfn.sparse, data = dsirData, times = times.d, pars = dde.fit$res$pars, beta = dde.fit$res$beta, active = NULL, coefs = dde.fit$res$coefs, basisvals = basis.d, lambda = 1000, in.meth='nlminb', basisvals0 = basis0, coefs0 = coefs0, nbeta = 2, ndelay = 2, tau = list(c(2, 7/3))))
    sim.Covar[[i]] <- Covar
    nnls.res[[i]] <- dde.fit$res
    save(nnls.res, sim.Covar, DSIR.pars, file =filename)
}

runTime <- Sys.time() - begTime
print(runTime)
save(nnls.res, sim.Covar, DSIR.pars, runTime, file = filename)
