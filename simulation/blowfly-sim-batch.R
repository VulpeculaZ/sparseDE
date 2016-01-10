library(gpDDE)
library(spam)
source("./R/make.blowfly.R")
## detach("package:limSolve", unload = TRUE)

load("data-blowfly-250-dist.RData")
blowfliesfn <- make.blowfly()

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
args <- commandArgs(TRUE)
dataRange <- (1 + 25 * as.numeric(args[1])) : (25 * (as.numeric(args[1]) + 1))
filename <- paste("blowfly-nnls-250-dist", as.numeric(args[1]),".RData", sep = "")
tau <- list(seq(5.5,10,0.5))


begTime <- Sys.time()
nnls.res <- list()
for(i in 1:length(dataRange)){
    print(dataRange[i])
    blowfly.data <- data.res[[dataRange[i]]]$blowfly.data
    blowfly.data.d <- blowfly.data[times0 >= 20]
    blowfly.data <- matrix(blowfly.data, length(blowfly.data),1)
    blowfly.data.d <- matrix(blowfly.data.d, length(blowfly.data.d),1)
    fdnames=list(NULL,c('y'),NULL)
    DEfd0 <- smooth.basis(times0, blowfly.data,fdPar(bbasis0,1,0.1))
    DEfd.d <- smooth.basis(times.d, blowfly.data.d,fdPar(bbasis.d,1,0.1))
    coefs0 <-  DEfd0$fd$coefs
    coefs.d <- DEfd.d$fd$coefs
    colnames(coefs0) <- colnames(coefs.d) <- c("y")
    ## plot the smooth plus data
    ## plotfit.fd(blowfly.data[,2],times0,DEfd0$fd)
    DEfd0 <- fd(coefs0,bbasis0, fdnames)
    DEfd.d <- fd(coefs.d,bbasis.d, fdnames)

    lambda <- 10000
    initBeta <- data.res[[dataRange[i]]]$initBeta
    initPars <- data.res[[dataRange[i]]]$initPars
    try(dde.fit <- Profile.LS.DDE(blowfliesfn, blowfly.data.d, times.d, pars = initPars, beta = initBeta, coefs = coefs.d, basisvals = bbasis.d, lambda = lambda, in.meth='nlminb',  basisvals0 = bbasis0, coefs0 = coefs0, nbeta = length(initBeta), ndelay = 1, tau = tau, control.out = list(method = "nnls.eq", maxIter = 50, lambda.sparse = 0, echo = TRUE)))
    nnls.res[[i]] <- dde.fit$res
    save(nnls.res, file =filename)
}

runTime <- Sys.time() - begTime
print(runTime)
save(nnls.res, runTime, file = filename)
