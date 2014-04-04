source("./R/timevar.R")
## source("./sources/LS.sparse.R")
## source("./sources/DSIRfnSparse.R")

library(penalized)
library(CollocInfer)
library(deSolve)
library(nnls)

## Function to simulate  SIR model:
tvSIR.gen <- function(t, y, parms){
    sint <- 8000 * (sin(t / parms["f"] / pi) / 2 + 2)
    dyS <- - (tvtrans(t, kappa)) * y[2] * y[1] +  sint
    dyI <- (tvtrans(t, kappa)) * y[2] * y[1] - parms["gamma"] * y[2]
    list(c(dyS, dyI))
}
yinit <- c(5000, 800)
times <- seq(0, 5, by = 1/52)
tvSIR.pars <- c(10, 1)
kappa <- rep(0.005, 12)
kappa[10:12] <- 0.002
names(tvSIR.pars) <- c("gamma", "f")
yout <- dede(y = yinit, times = times, func = tvSIR.gen, parms = tvSIR.pars, atol = 1e-10)
pdf("tvSIRsim.pdf", 9, 5)
matplot(yout[,1], yout[,-1], type = "l", lwd = 2, main = "Time Varying SIR Model", xlab="time", ylab="Numbers of individuals")
legend("topright", legend = c("S","I"), col=c(1,2), lty = c(1,2), lwd = c(2,2))
dev.off()
knots <- times
norder = 3
nbasis = length(knots) + norder - 2
range  = range(knots)
basis <- create.bspline.basis(range=range(knots), nbasis=nbasis, norder=norder, breaks=knots)
fdnames=list(NULL,c('S', 'I'),NULL)
bfdPar <- fdPar(basis,lambda=0.2,int2Lfd(1))
initUnif <- runif(100, -1,1)
initUnifKappa <- runif(100, -0.001,0.001)

begTime <- Sys.time()
set.seed(42)
sim.res <- list()
for(i in 1:100){
    print(i)
    xout <- c()
    xout <- cbind(xout, yout[,2] + rnorm(length(yout[,2]), sd = 100))
    xout <- cbind(xout, yout[,3] + rnorm(length(yout[,2]), sd = 100))
    ## points(times, xout)
    DEfd <- smooth.basis(knots, xout, bfdPar,fdnames=fdnames)$fd
    ## temp.fit <- eval.fd(times.d, DEfd.d)
    ## par(ask=FALSE)
    ## plotfit.fd(xout, times, DEfd)
    ## plotfit.fd(xout[times >=5,], times.d, DEfd.d)
    ## extract the coefficients and assign variable names
    coefs <- DEfd$coefs
    colnames(coefs) = c("S","I")
    ##  Set a value for lambda
    ## lambda = 1000
    ## Data
    tvData <- matrix(xout, ncol =2, dimnames = list(NULL,c("S", "I")))
    ## Setting initial values
    initPars <- 10 + initUnif[i]
    names(initPars) <- c("gamma")
    initKappa <- rep(0.005, 12)
    initKappa <- initKappa + initUnifKappa[i]
    names(initKappa) <- c("k1", "k2", "k3","k4","k5","k6","k7","k8","k9","k10","k11", "k12")
    tv.fit <- Profile.LS.tv(tvDSIRfn, tvData, times=times, pars = initPars, kappa = initKappa, coefs = coefs, basisvals = basis, lambda = 1000, in.meth='nlminb', control.out = list(method = "nnls", maxIter = 20, echo = FALSE, lambda.sparse = 0))
    sim.res[[i]] <- tv.fit
    curseed <- get(".Random.seed", .GlobalEnv)
    save(sim.res, curseed, file ="sim.tv02.sd100.RData")
}
curseed <- get(".Random.seed", .GlobalEnv)
runTime <- Sys.time() - begTime
print(runTime)
save(sim.res, runTime, curseed, file = "sim.tv02.sd100.RData")
