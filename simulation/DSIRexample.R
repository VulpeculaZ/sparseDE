library(deSolve)
library(CollocInfer)
library(MASS)
library(limSolve)
source("./R/poslasso.R")
source("./R/tv-delay-cov.R")
source("./R/sparse.R")
source("./R/LS.sparse.R")

library(gpDDE)
DSIR.gen <- function(t, y, parms){
    if(t < 0){
        lagI1 <- 2000
        lagI2 <- 2000
    }
    else{
        lagI1 <- lagvalue(t - parms["tau1"], 2)
        lagI2 <- lagvalue(t - parms["tau2"], 2)
    }
    sint <- 4000 * (sin(t / pi) + 2)
    dyS <- - parms["beta"] * (0.5 * ( lagI1 + lagI2)) * y[1] +  sint
    dyI <- parms["beta"] * (0.5 * ( lagI1 + lagI2)) * y[1] - parms["gamma"] * y[2]
    list(c(dyS, dyI))
}

set.seed(1234)
xinit <- c(4000, 2000)
DSIR.pars <- c(5, 0.0012, 0.6, 0.8)
names(DSIR.pars) <- c("gamma", "beta", "tau1", "tau2")
times <- seq(-DSIR.pars["tau2"], 30, by = 0.1)
xout <- dede(y = xinit, times = times, func = DSIR.gen, parms = DSIR.pars, atol = 1e-10)
yout <- c()
yout <- cbind(yout, xout[,2] + rnorm(length(xout[,2]), sd = 100))
yout <- cbind(yout, xout[,3] + rnorm(length(xout[,2]), sd = 100))
DSIRdata <- yout
save(DSIRdata, file = "DSIRdata.RData")
times <- seq(-DSIR.pars["tau2"], 30, by = 0.1)
yout0 <- yout[times >= 0, ]
yout.d <- yout[times >= 5, ]


pdf("DSIRsim.pdf",7, 4)
matplot(xout[,1], xout[,-1], type = "l", lwd = 2, main = "Delay SIR Model",  xlab = "Time", ylab="Population")
points(times0, yout0[ ,1])
points(times0, yout0[ ,2])
legend("topright", legend = c("S","I"), col=c(1,2), lty = c(2,2), lwd = c(2,2))
dev.off()


colnames(yout.d) <-  c("S","I")
times0 <- times[times>=0]
times.d <- times[times>=5]
norder = 3
nbasis.d = length(times.d) + norder - 2
nbasis0 <- length(times0) + norder - 2
basis0 <- create.bspline.basis(range=range(times0), nbasis=nbasis0, norder=norder, breaks=times0)
basis.d <- create.bspline.basis(range=range(times.d), nbasis=nbasis.d, norder=norder, breaks=times.d)
fdnames=list(NULL,c('S', 'I'),NULL)
bfdPar0 = fdPar(basis0,lambda=1,int2Lfd(1))
bfdPar.d <- fdPar(basis.d,lambda=1,int2Lfd(1))
DEfd0 <- smooth.basis(times0, yout0, bfdPar0,fdnames=fdnames)$fd
DEfd.d <- smooth.basis(times.d, yout.d, bfdPar.d, fdnames=fdnames)$fd
coefs0 <- DEfd0$coefs
colnames(coefs0) = c("S","I")
coefs.d <- DEfd.d$coefs
colnames(coefs.d) = c("S", "I")

initPars <- c(4, 0.001)
names(initPars) <- c("gamma", "beta")
initBeta <- rep(1/11, 11)
tau <- list(seq(0,2, length.out = 11))
##  Set a value for lambda
lambda = 1000
DSIRfn <- DSIRfn.make()


dde.fit <- Profile.LS.DDE(DSIRfn, yout.d, times.d, pars = initPars, beta = initBeta, coefs = coefs.d, basisvals = basis.d, lambda = 1000, in.meth='nlminb', basisvals0 = basis0, coefs0 = coefs0, nbeta = length(initBeta), ndelay = 2, tau = tau, control.out = list(method = "nnls.eq", maxIter = 50, lambda.sparse = 0, echo = TRUE))

DSIRInitCoefs <- dde.fit$ncoefs
save(DSIRInitCoefs, file = "DSIRinit.RData")

initPars <- c(5, 0.0012)
names(initPars) <- c("gamma", "beta")
initBeta <- rep(0, 11)
initBeta[c(4,5,11)] <- c(0.611, 0.362, 0.026)
tau <- list(seq(0,2, length.out = 11))
lambda = 1000
DSIRfn <- DSIRfn.make()
system.time(dde.fit <- Profile.LS.DDE(DSIRfn, yout.d, times.d, pars = initPars, beta = initBeta, coefs = DSIRInitCoefs, basisvals = basis.d, lambda = 1000, in.meth='nlminb', basisvals0 = basis0, coefs0 = coefs0, nbeta = length(initBeta), ndelay = 2, tau = tau, control.out = list(method = "nnls.eq", maxIter = 2, lambda.sparse = 0, echo = TRUE)))


pdf(file = "DSIRfit.pdf", width = 10, height = 5)
par(mfrow = c(1,2))
initFdobj.d <- fd(dde.fit$res$coefs, basis.d)
plotfit.fd(yout.d, times.d, DEfd.fit, xlab ="Time", ylab = "Population", ask = TRUE)
dev.off()

res <- yout.d - eval.fd(times.d, DEfd.fit)
res.arima <- arima(res[,2], include.mean = FALSE)
tsdiag.Arima(res.arima)

dde.fit1 <- dde.fit

DDEdiag(y = yout.d, times = times.d, fitted = DEfd.fit, use.TSA = TRUE)

IntegrateForward.DDE( times.forecast = c(25, 26), pars = dde.fit$res$pars, beta = dde.fit$res$beta,  proc = dde.fit$proc, more = dde.fit$more, tau = tau, ndelay = 2, fdobj0 = dde.fit$fdobj0, fdobj.d = DEfd.fit)

forecast.DDE(y = yout.d, times = times.d, h = 40, pars = dde.fit$res$pars, beta = dde.fit$res$beta, proc = dde.fit$proc, more = dde.fit$more, tau = tau, ndelay = 2, fdobj0 = dde.fit$fdobj0, fdobj.d = DEfd.fit, ask = FALSE, xlab = "times", ylab = "Population")

parms <- list(pars = dde.fit$res$pars, beta = dde.fit$res$beta,  proc = dde.fit$proc, more = dde.fit$more, tau = tau, ndelay = 2, fdobj0 = dde.fit$fdobj.d, fdobj.d = dde.fit$fdobj0, t0 = 25, tau.max = 0.6)
y0 <- c(eval.fd(25, dde.fit$fdobj.d))
out = dede(y = y0, times = c(25,26), func = dderhs, parms = parms)


sparse.fit <- sparse.DDE(fn = DSIRfn, data =  yout.d, times = times.d,
                        basisvals = basis.d, lambda = 1000, in.meth='nlminb',
                        basisvals0 = basis0, coefs0 = coefs0,
                        nbeta = 11, ndelay = 2, tau = tau,
                        control.out = list(method = "nnls", maxIter = 10,
                        selection.method = "lars"), nnls.res = dde.fit$res)

sparse.fit$select$beta

pdf(file = "DSIRlasso.pdf", width = 10, height = 5)
par(mfrow = c(1,2))
DEfd.fit <- fd(sparse.fit$select$coefs, basis.d)
plotfit.fd(yout.d, times.d, DEfd.fit, xlab = "Time",
           ylab = "Population")
dev.off()


Covar <- ProfileSSE.covariance.DDE(pars = dde.fit$res$pars,
                                     beta = dde.fit$res$beta, active = NULL,
                                     fn = DSIRfn, data = yout.d,
                                     times = times.d,  coefs = dde.fit$res$coefs,
                                     basisvals = basis.d, lambda = 1e7,
                                     in.meth='nlminb', delay = delay,
                                     basisvals0 = basis0, coefs0 = coefs0,
                                     nbeta = length(beta), ndelay = 2, tau = tau)

diag(Covar)

tvtrans <- function(t,k){
    period <- (t %% 5) / 5
    r <- rep(0, length(t))
    ka <- c(k, k[1])
    for(i in 1:length(k)){
        mk <- period[period  >=  (i-1)/length(k) & period <  i /length(k)]
        r[period  >=  (i-1)/length(k) & period <  i /length(k)] <-
            k[i] + (ka[i+1] - k[i]) * (mk - (i-1)/length(k)) * length(k)
    }
    return(r)
}


DTVSIR.gen <- function(t, y, parms, kappa){
    if(t < 0){
        lagI1 <- 2000
    }
    else{
        lagI1 <- lagvalue(t - parms["tau1"], 2)
    }
    sint <- 4000 * (sin(t / pi) + 2)
    dyS <- - (tvtrans(t, kappa)) * (lagI1) * y[1] +  sint
    dyI <- (tvtrans(t, kappa)) * (lagI1) * y[1] - parms["gamma"] * y[2]
    list(c(dyS, dyI))
}

xinit <- c(4000, 2000)
DTVSIR.pars <- c(5, 0.5)
kappa <- c(rep(0.002, 3), rep( 0.001, 3))

names(DTVSIR.pars) <- c("gamma", "tau1")
times <- seq(-DTVSIR.pars["tau1"], 30, by = 0.1)
xout <- dede(y = xinit, times = times, func = DTVSIR.gen, parms = DTVSIR.pars, atol = 1e-10, kappa = kappa)
yout <- c()
set.seed(42)

times <- seq(-0.5, 30, by = 0.1)
yout <- cbind(yout, xout[,2] + rnorm(length(xout[,2]), sd = 100))
yout <- cbind(yout, xout[,3] + rnorm(length(xout[,2]), sd = 100))
times0 <- times[times>=0]
times.d <- times[times>=5]
yout0 <- yout[times >= 0, ]
yout.d <- yout[times >= 5, ]


pdf("DTVSIRsim.pdf",7, 5)
matplot(xout[,1], xout[,-1], type = "l", lwd = 2, main = "Delay SIR Model",  xlab = "Time", ylab="Population")
points(times0, yout0[ ,1])
points(times0, yout0[ ,2])
legend("topright", legend = c("S","I"), col=c(1,2), lty = c(2,2), lwd = c(2,2))
dev.off()

colnames(yout.d) <-  c("S","I")
norder = 3
nbasis.d = length(times.d) + norder - 2
nbasis0 <- length(times0) + norder - 2
basis0 <- create.bspline.basis(range=range(times0), nbasis=nbasis0, norder=norder, breaks=times0)
basis.d <- create.bspline.basis(range=range(times.d), nbasis=nbasis.d, norder=norder, breaks=times.d)
fdnames=list(NULL,c('S', 'I'),NULL)
bfdPar0 = fdPar(basis0,lambda=1,int2Lfd(1))
bfdPar.d <- fdPar(basis.d,lambda=1,int2Lfd(1))
DEfd0 <- smooth.basis(times0, yout0, bfdPar0,fdnames=fdnames)$fd
DEfd.d <- smooth.basis(times.d, yout.d, bfdPar.d, fdnames=fdnames)$fd
coefs0 <- DEfd0$coefs
colnames(coefs0) = c("S","I")
coefs.d <- DEfd.d$coefs
colnames(coefs.d) = c("S", "I")

DTVSIRfn.make <- function(){
    fn <- function (t, y, p, more)
    {
        r = y
        yi.d <- more$y.d[,1]
        pk <- p[(length(p) - more$nKappa + 1):length(p)]
        r[, "S"] =  - tvtrans(t, pk) * yi.d * y[, "S"] + 4000 * (sin(t / pi) + 2)
        r[, "I"] =  tvtrans(t, pk) * yi.d * y[, "S"] - p["gamma"] * y[, "I"]
        return(r)
    }

    dfdx <- function (t, y, p, more)
    {
        r = array(0, c(length(t), ncol(y), ncol(y)))
        dimnames(r) = list(NULL, colnames(y), colnames(y))
        pk <- p[(length(p) - more$nKappa + 1):length(p)]
        yi.d <- more$y.d[,1]
        r[, "S", "S"] = - tvtrans(t, pk) * yi.d
        r[, "I", "S"] =  tvtrans(t, pk) * yi.d
        r[, "I", "I"] = -p["gamma"]
        return(r)
    }

    dfdx.d <- function (t, y, p, more)
    {
        pk <- p[(length(p) - more$nKappa + 1):length(p)]
        r = array(0, c(length(t), ncol(y), ncol(y)))
        dimnames(r) = list(NULL, colnames(y), colnames(y))
        r[, "S", "I"] = - tvtrans(t, pk) * y[,"S"]
        r[, "I", "I"] =  tvtrans(t, pk) * y[,"S"]
        return(r)
    }

    dfdp <- function (t, y, p, more)
    {
        yi.d <- more$y.d[,1]
        r = array(0, c(length(t), ncol(y), length(p)))
        dimnames(r) = list(NULL, colnames(y), names(p))
        r[, "I", "gamma"] = - y[, "I"]
        period <- t %% 1
        nKappa <- more$nKappa
        for(i in 1 : nKappa){
            r[ , "S", paste("k", i, sep ="")][period  >= (i-1)/nKappa & period < i /nKappa] =
                - (yi.d * y[, "S"])[period  >= (i-1)/nKappa & period < i /nKappa]
            r[ , "I", paste("k", i, sep ="")][period  >= (i-1)/nKappa & period < i /nKappa] =
                (yi.d * y[, "S"])[period  >= (i-1)/nKappa & period < i /nKappa]
        }
        return(r)
    }

    d2fdx2 <- function (t, y, p, more)
    {
        r = array(0, c(length(t), ncol(y), ncol(y), ncol(y)))
        dimnames(r) = list(NULL, colnames(y), colnames(y), colnames(y))
        return(r)
    }

    d2fdxdp <- function (t, y, p, more)
    {
        yi.d <- more$y.d[,1]
        r = array(0, c(length(t), ncol(y), ncol(y), length(p)))
        dimnames(r) = list(NULL, colnames(y), colnames(y), names(p))
        r[, "I", "I", "gamma"] = -1
        period <- (t %% 5) / 5
        nKappa <- more$nKappa
        for(i in 1:nKappa){
            r[ , "S", "S", paste("k", i, sep ="")][period  >= (i-1)/nKappa & period < i /nKappa] = - yi.d[period  >= (i-1)/nKappa & period < i /nKappa]
            r[ , "I", "S", paste("k", i, sep ="")][period  >= (i-1)/nKappa & period < i /nKappa] = yi.d[period  >= (i-1)/nKappa & period < i /nKappa]
        }
        return(r)
    }

    d2fdx.ddp <- function (t, y, p, more)
    {
        yi.d <- more$y.d[,1]
        r = array(0, c(length(t), ncol(y), ncol(y), length(p)))
        dimnames(r) = list(NULL, colnames(y), colnames(y), names(p))
        period <- (t %% 5) / 5
        nKappa <- more$nKappa
        for(i in 1:nKappa){
            r[ , "S", "I", paste("k", i, sep ="")][period  >= (i-1)/nKappa & period < i /nKappa] = - y[,"S"][period  >= (i-1)/nKappa & period < i /nKappa]
            r[ , "I", "I", paste("k", i, sep ="")][period  >= (i-1)/nKappa & period < i /nKappa] = y[,"S"][period  >= (i-1)/nKappa & period < i /nKappa]
        }
        return(r)
    }


    d2fdxdx.d <- function (t, y, p, more)
    {
        yi.d <- more$y.d[,1]
        pk <- p[(length(p) - more$nKappa + 1):length(p)]
        r = array(0, c(length(t), ncol(y), ncol(y), ncol(y)))
        dimnames(r) = list(NULL, colnames(y), colnames(y), colnames(y))
        r[, "S", "S", "I"] = - tvtrans(t, pk)
        r[, "I", "S", "I"] =  tvtrans(t, pk)
        return(r)
    }


    d2fdx.d2 <- function (t, y, p, more)
    {
        yi.d <- more$y.d[,1]
        r = array(0, c(length(t), ncol(y), ncol(y), ncol(y)))
        dimnames(r) = list(NULL, colnames(y), colnames(y), colnames(y))
        return(r)
    }

    return(list(
        fn = fn, dfdx = dfdx,
        dfdp = dfdp, d2fdx2 = d2fdx2,
        d2fdxdp = d2fdxdp, d2fdx.ddp = d2fdx.ddp,
        dfdx.d = dfdx.d, d2fdx.ddx = d2fdxdx.d,
        d2fdxdx.d = d2fdxdx.d, d2fdx.d2 = d2fdx.d2
        ))
}

DTVSIRfn <- DTVSIRfn.make()


initPars <- c(4)
names(initPars) <- c("gamma")
initBeta <- rep(0.2,5)
initKappa <- c(rep(0.0015, 3), rep(0.0015,3))
names(initKappa) <- c("k1", "k2", "k3","k4" ,"k5","k6")
tau <- list(seq(0,2, length.out = 5))
##  Set a value for lambda
lambda = 1000

DTVSIRfn <- DTVSIRfn.make()
source("./R/tv-delay.R")

dde.fit <- Profile.LS.tv.delay(DTVSIRfn, yout.d, times.d, pars = initPars, beta = initBeta, kappa = initKappa, coefs = coefs.d, basisvals = basis.d, lambda = 1000, in.meth='nlminb', delay = delay, basisvals0 = basis0, coefs0 = coefs0, nbeta = length(initBeta), ndelay = 2, tau = tau, control.out = list(method = "nnls.eq", maxIter = 50, lambda.sparse = 0, echo = TRUE))
