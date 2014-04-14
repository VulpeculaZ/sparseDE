library(CollocInfer)
source("./R/SIRfn.R")
source("./R/dsDE.R")
## Function to simulate  SIR model:
SIR.gen <- function(t, y, parms){
    sint <- sin(t / parms["f"]) + 1
    dyS <- -parms["beta"] * y[2] * y[1] + (1 - y[1]) * (parms["b"]* sint + 0.1)
    dyI <- parms["beta"] * y[2] * y[1] - (parms["gamma"] + (parms["b"] * sint + 0.1)) * y[2]
    list(c(dyS, dyI))
}



yinit <- c(0.95, 0.05)
times <- seq(0, 25, by = 0.1)
## Simulation for SIR
SIR.pars <- c(2, 0.25, 0.5, 2)
names(SIR.pars) <- c("beta", "b","gamma", "f")
yout <- ode(y = yinit, times = times, func = SIR.gen, parms = SIR.pars, atol = 1e-10)
matplot(yout[,1], yout[,-1], type = "l", lwd = 2, main = "SIR Model")


knots <- times
norder = 3
nbasis = length(knots) + norder - 2
range  = range(knots)
basis <- create.bspline.basis(range=range(knots), nbasis=nbasis, norder=norder, breaks=knots)
fdnames=list(NULL,c('S', 'I'),NULL)
bfdPar = fdPar(basis,lambda=1,int2Lfd(1))

set.seed(42)
xout <- c()
xout <- cbind(xout, yout[,2] + rnorm(length(yout[,2]), sd = 0.01))
xout <- cbind(xout, yout[,3] + rnorm(length(yout[,2]), sd = 0.01))
## points(times, xout)
DEfd <- smooth.basis(knots, xout, bfdPar,fdnames=fdnames)$fd
## temp.fit <- eval.fd(times.d, DEfd.d)
## par(ask=FALSE)
plotfit.fd(xout, times, DEfd)

## extract the coefficients and assign variable names
coefs <- DEfd$coefs
colnames(coefs) = c("S","I")
##  Set a value for lambda
lambda = 1000
## Data
dsirData <- matrix(xout, ncol =2, dimnames = list(NULL,c("S", "I")))

## Setting initial values
initPars <- c(0.6, 2)
names(initPars) <- c("gamma", "beta")



## Estimation:

## One beta
debug(Profile.LS)
system.time(dde.fit <- Profile.LS(SIRfn, dsirData, times, initPars, coefs = coefs, basisvals = basis, lambda, in.meth='nlminb'))
=DEfd2 = fd(dde.fit$coefs,basis.d, fdnames)
plotfit.fd(dsirData, times.d , DEfd2)


debug(dsde.LS)
system.time(dsDE.fit <- dsde.LS(SIRfn, dsirData, times, initPars, coefs = coefs, basisvals = basis, lambda = 1000))

## 3 betas
initBeta <- c(1,0.5,0.5)
dde.fit <- Profile.LS.sparse(DSIRfn.sparse, dsirData, times.d, initPars,beta = initBeta, coefs = coefs.d, basisvals = basis.d, lambda, in.meth='nlminb', delay = delay, basisvals0 = basis0, coefs0 = coefs0, nbeta = 3, ndelay = 2, tau = list(c(0,1,2)), control.out = list(method = "nnls", maxIter = 20, lambda.sparse = 0))
DEfd2 = fd(dde.fit$coefs,basis.d, fdnames)
plotfit.fd(dsirData, times.d , DEfd2)

source("./sources/poslasso.R")
## start with nnls converged values:
debug(nls.sparse)
dde.fit <- Profile.LS.sparse(DSIRfn.sparse, dsirData, times.d, pars = dde.fit$pars, beta = dde.fit$beta, coefs = dde.fit$coefs, basisvals = basis.d, lambda = 100, in.meth='nlminb', delay = delay, basisvals0 = basis0, coefs0 = coefs0, nbeta = 3, ndelay = 2, tau = list(c(0,1,2)), control.out = list(method = "nnls", maxIter = 20, lambda.sparse = 0.1, tol = 1e-6))
tmp <- lars.pos(x= Zdf, y= y)

## 5 betas
initBeta <- c(0.2,1.2,0.2, 0.2, 0.2)
dde.fit <- Profile.LS.sparse(DSIRfn.sparse, dsirData, times.d, initPars,beta = initBeta, coefs = coefs.d, basisvals = basis.d, lambda, in.meth='nlminb', delay = delay, basisvals0 = basis0, coefs0 = coefs0, nbeta = 5, ndelay = 2, tau = list(c(0,1,2,3,4)), control.out = list(method = "nnls", maxIter = 20, lambda.sparse = 0))
## start with nnls converged values:
dde.fit <- Profile.LS.sparse(DSIRfn.sparse, dsirData, times.d, pars = dde.fit$pars, beta = dde.fit$beta, coefs = dde.fit$coefs, basisvals = basis.d, lambda = 100, in.meth='nlminb', delay = delay, basisvals0 = basis0, coefs0 = coefs0, nbeta = 5, ndelay = 2, tau = list(c(0,1,2,3,4)), control.out = list(method = "nnls", maxIter = 1, lambda.sparse = 0))

## 9 betas:
initBeta <- c()
for(i in 1:(length(dde.fit$beta)-1)){
    initBeta <- c(initBeta, dde.fit$beta[i], (dde.fit$beta[i] + dde.fit$beta[i+1])/2)
}
initBeta <- c(initBeta, dde.fit$beta[5])
initBeta <- initBeta * 5 / 9
## dde.fit <- Profile.LS.sparse(DSIRfn.sparse, dsirData, times.d, pars = dde.fit$pars, beta = initBeta, coefs = dde.fit$coefs, basisvals = basis.d, lambda = 100, in.meth='nlminb', delay = delay, basisvals0 = basis0, coefs0 = coefs0, nbeta = length(initBeta), ndelay = 2, tau = list(seq(0,4, length.out =9)), control.out = list(method = "nnls", maxIter = 20, lambda.sparse = 0))

lambda0 = max(abs(as.vector(t(dde.fit$outer.result$y) %*% dde.fit$outer.result$Zdf)))
lambda = exp(seq(log(lambda0), log(lambda0 * 0.001), len = 20))

lambda.min <- lambda[10]

b9.fit <- list()
for(i in 1:length(lambda)){
    b9.fit[[i]] <- Profile.LS.sparse(DSIRfn.sparse, dsirData, times.d, pars = dde.fit$pars, beta = initBeta, coefs = dde.fit$coefs, basisvals = basis.d, lambda = 100, in.meth='nlminb', delay = delay, basisvals0 = basis0, coefs0 = coefs0, nbeta = length(initBeta), ndelay = 2, tau = list(seq(0,4, length.out =9)), control.out = list(method = "penalized", maxIter = 20, lambda.sparse = lambda.min[i]))
}

## 3 betas, using penalized
initBeta <- c(1,0.5,0.5)
initPars <- c(0.500)
names(initPars) <- c("gamma")
## debugonce(nls.sparse)
dde.fit <- Profile.LS.sparse(DSIRfn.sparse, dsirData, times.d, initPars,beta = initBeta, coefs = coefs.d, basisvals = basis.d, lambda = 100, in.meth='nlminb', delay = delay, basisvals0 = basis0, coefs0 = coefs0, nbeta = 3, ndelay = 2, tau = list(c(0,1,2)), control.out = list(method = "penalized", maxIter = 20, lambda.sparse = 0.1, tol = 1e-6))


DEfd2 = fd(dde.fit$coefs,basis.d, fdnames)
plotfit.fd(dsirData, times.d , DEfd2)


## Uses two-stage:
debugonce(nls.sparse)
initBeta <- c(0.6,0.6,0.6)
initPars <- c(0.600)
names(initPars) <- c("gamma")
dde.fit <- Profile.LS.sparse(DSIRfn.sparse, dsirData, times.d, initPars,beta = initBeta, coefs = coefs.d, basisvals = basis.d, lambda = 100, in.meth='nlminb', delay = delay, basisvals0 = basis0, coefs0 = coefs0, nbeta = 3, ndelay = 2, tau = list(c(0,1,2)), control.out = list(method = "twoStage", maxIter = 20, lambda.sparse = 0.3, tol = 1e-6))
DEfd2 = fd(dde.fit$coefs,basis.d, fdnames)
plotfit.fd(dsirData, times.d , DEfd2)

source("./sources/SIRfn.R")
lambda = 1000
SIRinitPars <- c(0.6, 2)
names(SIRinitPars) <- c("gamma", "beta")
debugonce(Profile.LS)
sir.fit <- Profile.LS(SIRfn, dsirData, times.d, SIRinitPars, coefs = coefs.d, basisvals = basis.d, lambda, in.meth='nlminb')

allpars = pars
allpars[active] = pars
coefs <- ncoefs
f <- ProfileSSE.AllPar(pars = allpars, times = times, data = data, coefs = coefs, lik = lik, proc = proc, in.meth = in.meth, control.in = control.in, dcdp = NULL, oldpars = NULL, use.nls = TRUE, sgn = 1)
## Check orcal and compare it with it.

Zdf <- res$Zdf
Xdf <- res$Xdf
y <- res$y
lambda0 = max(abs(as.vector(t(y) %*% Zdf)))
lambda = exp(seq(log(lambda0), log(lambda0 * 0.001), len = 100))
delta <- rep(1, length(y))
pars.ts.mcp <- beta.ts.mcp <- pars.ts.lasso <- beta.ts.lasso <- pars.pen <- beta.pen <- c()

for(i in 1:length(lambda)){
    lambda.sparse <- lambda[i]
    ## res <- twoStgEst(Z = Zdf, y = y, X = Xdf, delta = delta, n_lambda = 1, lambda0 = lambda.sparse)
    ## pars.ts.mcp <- rbind(pars.ts.mcp,res$betamcp)
    ## beta.ts.mcp <- rbind(beta.ts.mcp, res$thetamcp)
    ## pars.ts.lasso <- rbind(pars.ts.lasso, res$betalasso)
    ## beta.ts.lasso <- rbind(beta.ts.lasso, res$thetalasso)
    res.sparse <- penalized(response = y, penalized = Zdf, unpenalized = Xdf, lambda1 = lambda.sparse)
    pars.pen <- rbind(pars.pen, res.sparse@unpenalized)
    beta.pen <- rbind(beta.pen, res.sparse@penalized)
}

#res <- lm.fit(x=cbind(Xdf, Zdf), y = y)
res <- lars.pos(x= cbind(Xdf,Zdf), y= y, positive = TRUE)
res$beta
plot(res)

