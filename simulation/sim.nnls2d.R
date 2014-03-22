source("./sources/sparse.R")
source("./sources/LS.sparse.R")
source("./sources/DSIRfnSparse.R")
source("./sources/poslasso.R")

library(penalized)
library(CollocInfer)
library(deSolve)
library(nnls)

## Function to simulate DSIR
DSIR.gen <- function(t, y, parms){
    if(t < 0){
        lagI1 <- 0.05
        lagI2 <- 0.05
    }
    else{
        lagI1 <- lagvalue(t - parms["tau1"], 2)
        lagI2 <- lagvalue(t - parms["tau2"], 2)
    }
    sint <- sin(t / parms["f"]) + 1
    if(y[2] > 1)
        y[2] <- 1
    if(y[2] < 0)
        y[2] <- 0
    dyS <- -parms["beta"] * (lagI1 + lagI2) * y[1] + (1 - y[1]) * (parms["b"] * sint + 0.1)
    dyI <- parms["beta"] * (lagI1 + lagI2) * y[1] - (parms["gamma"] + (parms["b"] * sint + 0.1)) * y[2]
    list(c(dyS, dyI))
}

yinit <- c(0.9, 0.1)
times <- seq(0, 25, by = 0.1)
SIR.pars <- c(2, 0.25, 0.5, 1)
names(SIR.pars) <- c("beta", "b","gamma", "f")
## Simulation for DSIR
tau <- c(1,2)
DSIR.pars <- c(SIR.pars, tau)
names(DSIR.pars) <- c("beta", "b","gamma", "f", "tau1", "tau2")
times <- seq(-DSIR.pars["tau2"], max(times), by = 0.1)
yout <- dede(y = yinit, times = times, func = DSIR.gen, parms = DSIR.pars, atol = 1e-10)

pdf("DSIRsim.pdf",7, 4)
matplot(yout[21:271,1], yout[21:271,-1], type = "l", lwd = 2, main = "Delayed SIR Model", xlab="time", ylab="Proportion of Population")
legend("topright", legend = c("S","I"), col=c(1,2), lty = c(1,2), lwd = c(2,2))
dev.off()

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

set.seed(42)
sim.res <- list()
for(i in 1:100){

    xout <- c()
    xout <- cbind(xout, yout[,2] + rnorm(length(yout[,2]), sd = 0.01))
    xout <- cbind(xout, yout[,3] + rnorm(length(yout[,2]), sd = 0.01))
    ## points(times, xout)
    xout0 <- xout[times >= 0,]
    xout.d <- xout[times >= 5,]
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
    initPars <- c(0.55)
    names(initPars) <- c("gamma")
    initBeta <- rep(0.1, 16)
    initBeta[4] <- initBeta[7] <- 0.6
    dde.fit <- Profile.LS.sparse(DSIRfn.sparse, dsirData, times.d, pars = initPars, beta = initBeta, coefs = coefs.d, basisvals = basis.d, lambda = 100, in.meth='nlminb', delay = delay, basisvals0 = basis0, coefs0 = coefs0, nbeta = length(initBeta), ndelay = 2, tau = list(seq(0,5, length.out = 16)), control.out = list(method = "nnls", maxIter = 10, lambda.sparse = 0))
    sim.res[[i]] <- dde.fit
    save(sim.res, DSIR.pars, file ="nnls2d01.RData")

}

pdf("DSIRsim.pdf",7, 4)
plot(times.d, dsirData[,1], ylim = c(0,0.5), col =1, xlab = "Time", ylab="Proportion of Population")
points(times.d, dsirData[,2], col=2)
tmp <- eval.fd(seq(5,25,length.out = 400), DEfd2)
lines(seq(5,25,length.out = 400), tmp[,1], lwd = 2, lty =2, col= 1)
lines(seq(5,25,length.out = 400), tmp[,2], lwd = 2, lty =2, col= 2)
legend("topright", legend = c("S","I"), col=c(1,2), lty = c(2,2), lwd = c(2,2))
dev.off()

DEfd2 = fd(dde.fit$ncoefs, basis.d, fdnames)
plotfit.fd(dsirData, times.d , DEfd2)
