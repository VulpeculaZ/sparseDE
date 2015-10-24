library(CollocInfer)
library(deSolve)
source("./R/vector.D.funcs.R")
source("./R/sparse.R")
source("./R/LS.sparse.R")
source("./R/poslasso.R")
source("./R/tv-delay-cov.R")

library(limSolve)
library(trustOptim)
library(MASS)
## detach("package:limSolve", unload = TRUE)


## Function to simulating data
vectorD.Gen <- function(t, y, parms){
    if(t<0)
        lag <- 0.05
    else
        lag <- lagvalue(t - parms["tau"])
    dy <- (parms["b"] + sin(t)) * lag * (1 - y) - parms["a"]* y
    list(dy, dy = dy)
}

vectorPars <- c(1, 1.8,0.8)
names(vectorPars) <- c("a","b","tau")
yinit <- 0.05
times <- seq(-0.8, 25, by = 0.1)
yout <- dede(y = yinit, times = times, func = vectorD.Gen, parms = vectorPars, atol = 1e-10)


times0 <- knots0 <- times[times >= 0]
times.d <- knots.d <- times[times >= 5]
norder = 3
nbasis.d = length(knots.d) + norder - 2
nbasis0 <- length(knots0) + norder - 2
range0  = range(knots0)
range.d <- range(knots.d)

## vectorPars <- c(1,1.8,0.8)
## names(vectorPars) <- c("a","b","tau")
vectorBasis0 <- create.bspline.basis(range=range(knots0), nbasis=nbasis0, norder=norder, breaks=knots0)
vectorBasis.d <- create.bspline.basis(range=range.d, nbasis=nbasis.d, norder=norder, breaks=knots.d)

##  Set a value for lambda
lambda = 1000
fdnames=list(NULL,c('S'),NULL)
bfdPar0 = fdPar(vectorBasis0,lambda=1,int2Lfd(1))
bfdPar.d <- fdPar(vectorBasis.d,lambda=1,int2Lfd(1))

## Setting initial values
vectorPars <- c(0.8, 1.5)
names(vectorPars) <- c("a","b")
tau <- list(seq(0, 2, by = 0.2))
beta <- rep(0, 11)
beta[4:8] <- 0.2

xout <- yout[,2] + rnorm(length(yout[,2]), sd = 0.01)

xout.d <- xout[times >= 5]
xout[xout < 0] <- 0
xout0 <- xout[times >= 0]
DEfd0 <- smooth.basis(knots0,xout0,bfdPar0,fdnames=fdnames)
DEfd.d <- smooth.basis(knots.d, xout.d, bfdPar.d, fdnames=fdnames)
## temp.fit <- eval.fd(times.d, DEfd.d)
## par(ask=FALSE)
## plotfit.fd(xout[times >=0], times0, DEfd0)
## plotfit.fd(xout[times >=5], times.d, DEfd.d)
## extract the coefficients and assign variable names
coefs0 <- DEfd0$fd$coefs
colnames(coefs0) = "S"
coefs.d <- DEfd.d$fd$coefs
colnames(coefs.d) = "S"
## Data
vectorData <- matrix(xout.d, ncol =1, dimnames = list(NULL,"S"))

vectorFun <- make.vector.disease.fn()


system.time(dde.fit <- Profile.LS.sparse(fn = vectorFun, data = vectorData, times = times.d, pars = vectorPars,  beta = beta, coefs = coefs.d, basisvals = vectorBasis.d,  lambda = 1e7, in.meth='nlminb', delay = delay, basisvals0 = vectorBasis0, coefs0 = coefs0, nbeta = length(beta), ndelay = 1, tau = tau, control.out = list(method = "nnls", maxIter = 20, lambda.sparse = 0, echo = TRUE, tol = 1e-12)))

sparse.fit <- LS.sparse(fn = vectorFun, data =  vectorData, times = times.d,
                        basisvals = vectorBasis.d, lambda = 1e7, in.meth='nlminb',
                        delay = delay, basisvals0 = vectorBasis0, coefs0 = coefs0,
                        nbeta = length(beta), ndelay = 1, tau = tau,
                        control.out = list(method = "nnls", maxIter = 10,
                        lambda.sparse = -3), nnls.res = dde.fit$res)

sparse.fit$select$beta

pdf(file = "vector-sim-lasso.pdf", width = 7, height = 5)
DEfd.fit <- fd(sparse.fit$select$coefs, vectorBasis.d)
plotfit.fd(vectorData, times.d, DEfd.fit, xlab = "Time",
           ylab = "Infection",  main = "Vector Disease Model")
dev.off()

Covar <- ProfileSSE.covariance.delay(pars = dde.fit$res$pars,
                                     beta = dde.fit$res$beta, active = NULL,
                                     fn = vectorFun, data = vectorData,
                                     times = times.d,  coefs = dde.fit$res$coefs,
                                     basisvals = vectorBasis.d, lambda = 1e7,
                                     in.meth='nlminb', delay = delay,
                                     basisvals0 = vectorBasis0, coefs0 = coefs0,
                                     nbeta = length(beta), ndelay = 1, tau = tau)
diag(Covar)
