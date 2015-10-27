begTime <- Sys.time()
source("./R/sparse.R")
source("./R/LS.sparse.R")
source("./R/DSIRfnSparse.R")
source("./R/poslasso.R")
source("./R/tv-delay-cov.R")
library(CollocInfer)
library(MASS)
library(spam)
library(nnls)
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
nnlsname <- paste("nnls-2dadj-sd02-", as.numeric(args[1]),".RData", sep = "")
filename <- paste("pen-2dadj-sd02-", as.numeric(args[1]),".RData", sep = "")
load(nnlsname)

set.seed(42)
sim.res.lars <- list()
sim.res.adlars <- list()
sim.Covar <- list()
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
    dde.fit <- LS.sparse(DSIRfn.sparse, data = dsirData ,times.d, basisvals = basis.d, lambda = lambda, in.meth='nlminb', delay = delay, basisvals0 = basis0, coefs0 = coefs0, nbeta = 16, ndelay = 2, tau = list(seq(0,5, length.out = 16)), control.out = list(lambda.sparse = -3), nnls.res = nnls.res[[i]])
    dde.fit2 <- LS.sparse(DSIRfn.sparse, data = dsirData ,times.d, basisvals = basis.d, lambda = lambda, in.meth='nlminb', delay = delay, basisvals0 = basis0, coefs0 = coefs0, nbeta = 16, ndelay = 2, tau = list(seq(0,5, length.out = 16)), control.out = list(lambda.sparse = -4), nnls.res = nnls.res[[i]])
    Covar <- ProfileSSE.covariance.delay(fn = DSIRfn.sparse, data = dsirData, times = times.d, pars = nnls.res[[i]]$pars, beta = nnls.res[[i]]$beta, active = NULL, coefs = nnls.res[[i]]$coefs, basisvals = basis.d, lambda = 1000, in.meth='nlminb', delay = delay, basisvals0 = basis0, coefs0 = coefs0, nbeta = length(nnls.res[[i]]$beta), ndelay = 2, tau = list(seq(0,5, length.out = 16)))
    sim.res.adlars[[i]] <- dde.fit$select
    sim.res.lars[[i]] <- dde.fit2$select
    sim.Covar[[i]] <- Covar
    save(sim.res.adlars, sim.res.lars, sim.Covar, file = filename)
}

runTime <- Sys.time() - begTime
scriptName <- "pen-16d-2dadj.R"
print(runTime)
save(sim.res.adlars, sim.res.lars, sim.Covar, DSIR.pars, runTime, scriptName, file = filename)
