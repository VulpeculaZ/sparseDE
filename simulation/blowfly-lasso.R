library(gpDDE)
library(spam)

args <- commandArgs(TRUE)
dataRange <- (1 + 25 * as.numeric(args[1])) : (25 * (as.numeric(args[1]) + 1))
nnls.filename <- paste("blowfly-nnls-250-", as.numeric(args[1]),".RData", sep = "")
res.filename <- paste("blowfly-lasso-250-", as.numeric(args[1]),".RData", sep = "")
dataRange <- (1 + 25 * as.numeric(args[1])) : (25 * (as.numeric(args[1]) + 1))

load(nnls.filename)
load("data-blowfly-250.RData")


blowfly.day <- seq(0,175, 0.5)
rr     = range(blowfly.day)       #  the range of observations times
knots  = seq(rr[1],rr[2],0.5)  #  knots at equally spaced values
norder = 3                      #  the order of the B-spline basis functions,
                                #  in this case piece-wise quadratic
nbasis = length(knots)+norder-2 #  the number of basis functions

#  set up the basis object
bbasis0 <- create.bspline.basis(range=rr, norder=norder, nbasis=nbasis,breaks=knots)
times0  <- blowfly.day
times.d  <- blowfly.day[blowfly.day >= 20]
knots.d <- seq(20,rr[2],0.5)
nbasis.d <- length(knots.d) + norder - 2
bbasis.d <- create.bspline.basis(range=c(20,rr[2]), norder=norder, nbasis=nbasis.d, breaks=knots.d)
bfdPar0 = fdPar(bbasis0,lambda=1,int2Lfd(1))
bfdPar.d <- fdPar(bbasis.d,lambda=1,int2Lfd(1))
tau <- list(seq(5.5,10,0.5))


## xout0 <- data.res[[1]]$xout[times >= 0.5,]
## DEfd0 <- smooth.basis(knots0, xout0, bfdPar0,fdnames=fdnames)$fd
## coefs0 <- DEfd0$coefs

set.seed(42)
sim.res.lars <- list()
sim.res.adlars <- list()

for(i in 1:length(nnls.res)){
    print(dataRange[i])
    initPars <- nnls.res[[i]]$pars
    initBeta <- nnls.res[[i]]$beta
    blowfly.data <- data.res[[dataRange[i]]]$blowfly.data
    blowfly.data.d <- blowfly.data[times0 >= 20]
    blowfly.data <- matrix(blowfly.data, length(blowfly.data),1)
    blowfly.data.d <- matrix(blowfly.data.d, length(blowfly.data.d),1)
    fdnames=list(NULL,c('y'),NULL)
    DEfd0 <- smooth.basis(times0, blowfly.data,fdPar(bbasis0,1,0.1))
    coefs0 <-  DEfd0$fd$coefs
    dde.fit <- sparse.DDE(blowfliesfn, data =  blowfly.data.d, times.d, basisvals = bbasis.d, lambda = 1000, in.meth='nlminb', basisvals0 = bbasis0, coefs0 = coefs0, nbeta = length(initBeta), ndelay = 1, tau = tau, control.out = list(method = "nnls.eq", maxIter = 20, selection.method = "lars"), nnls.res = nnls.res[[i]])
    sim.res.adlars[[i]] <- dde.fit$select
    dde.fit <- sparse.DDE(blowfliesfn, data =  blowfly.data.d, times.d, basisvals = bbasis.d, lambda = 1000, in.meth='nlminb', basisvals0 = bbasis0, coefs0 = coefs0, nbeta = length(initBeta), ndelay = 1, tau = tau, control.out = list(method = "nnls.eq", maxIter = 20, selection.method = "addaptive"), nnls.res = nnls.res[[i]])
    sim.res.lars[[i]] <- dde.fit$select
    save(sim.res.lars, sim.res.adlars, file = res.filename)
}


