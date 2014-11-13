source("./R/vector.D.funcs.R")
library(deSolve)
source("./R/sparse.R")
library(limSolve)

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
## solve the model

yout <- dede(y = yinit, times = times, func = vectorD.Gen, parms = vectorPars, atol = 1e-10)
# plot(yout, which = 1, type = "l", lwd = 2, main = "Vector Disease Model")
# plot(yout[,2], yout[,3], xlab = "y", ylab = "dy", type = "l", lwd = 2)


vect.time <- times[times >= 5]
rr     = range(vect.time)       #  the range of observations times
knots  = seq(rr[1],rr[2],2)  #  knots at 52 equally spaced values
norder = 3                      #  the order of the B-spline basis functions,
                                #  in this case piece-wise quadratic
nbasis = length(knots)+norder-2 #  the number of basis functions

#  set up the basis object

bbasis0 <- create.bspline.basis(range=rr, norder=norder, nbasis=nbasis,breaks=knots)
times0  <- vect.time
times.d  <- vect.time[vect.time >= 8]
knots.d <- seq(times.d[1] ,rr[2],0.1)
nbasis.d <- length(knots.d) + norder - 2
bbasis.d <- create.bspline.basis(range=c(times.d[1],rr[2]), norder=norder, nbasis=nbasis.d, breaks=knots.d)


## Generating Data
vec.data <- matrix(yout[,2][yout[,1] >=5 ] + rnorm(sum(yout[,1] >=5), sd = 0.01), length(times0),1)
colnames(vec.data) <- c("S")
vec.data.d <- vec.data[times0 >= 8, , drop=FALSE]

# To get an initial estimate of the states we smooth the observed I component
# and set the other coefficients to zero.

# smooth the log I values
fdnames=list(NULL,c('S'),NULL)
DEfd0 <- smooth.basis(times0 ,(vec.data),fdPar(bbasis0,1,0.1))
DEfd.d <- smooth.basis(times.d, vec.data.d,fdPar(bbasis.d,1,0.1))
coefs0 <-  DEfd0$fd$coefs
coefs.d <- DEfd.d$fd$coefs
colnames(coefs0) <- colnames(coefs.d) <- c("S")
# set up the functional data object for the three variables
# plot the smooth plus data
## plotfit.fd(vec.data[,2],times0,DEfd0$fd)
DEfd0 <- fd(coefs0,bbasis0, fdnames)
DEfd.d <- fd(coefs.d,bbasis.d, fdnames)


## Setting initial values
vectorPars <- c(0.8, 1.5)
names(vectorPars) <- c("a","b")
initBeta <- rep(0, 10)
initBeta[5:7] <- 1/5
initBeta[6] <- 3/5


dde.fit <- Profile.LS.sparse(vectorFun, vec.data.d, times.d, pars = vectorPars, beta = initBeta, coefs = coefs.d, basisvals = bbasis.d,  lambda=1000, in.meth='nlminb', delay = delay, basisvals0 = bbasis0, coefs0 = coefs0,  nbeta = length(initBeta), ndelay = 1, tau = list(seq(0.3,1.2, by = 0.1)), control.out = list(method = "nnls", maxIter = 20, lambda.sparse = 0))

DEfd2 = fd(dde.fit$coefs,vectorBasis.d, fdnames)
plotfit.fd(vectorData,times.d , DEfd2)
