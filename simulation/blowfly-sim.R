library(deSolve)
source("./R/blowflies.R")
source("./R/sparse.R")
library(limSolve)
## detach("package:limSolve", unload = TRUE)

init.beta.zoom <- function(beta){
    beta.new <- rep(NA, 2 * length(beta) -1 )
    for(i in 1:length(beta.new)){
        beta.new[i] <- ifelse(i %% 2 == 0, 0.25 * beta[i/2]+ 0.25 * beta[i/2+1], 0.5 * beta[(i+1)/2])
    }
    beta.new <- beta.new / sum(beta.new)
    beta.new
}

blowfly.gen <- function(t, y, parms){
    if(t<0)
        lag <- 1000
    else
        lag <- lagvalue(t - parms["tau"])
    dy <- parms["c"] * lag * exp(-lag / parms["N0"]) - parms["a"] * y
    list(dy, dy = dy)
}

blowfly.pars <- c(150 / 8, 8 / 8 , 1000, 15)
names(blowfly.pars) <- c("c", "a", "N0", "tau")
times <- seq(-10, 200, by = 0.5)

yout <- dede(y = 1000, times = times, func = blowfly.gen, parms = blowfly.pars, atol = 1e-7)
plot(yout[55:405,], which = 1, type = "l", lwd = 2, main = "Nicholson's Blowflies Model")
blowfly.day <- yout[,1][55:405] - yout[55,1]
blowfly.data <- yout[,2][55:405]   + rnorm(351, sd = 500)

pdf(file = "blowfly-sim.pdf", width = 7, height = 5)
plot(x = blowfly.day, y = yout[,2][55:405], type= "l", ylim = c(-500,8000), main = "Nicholson's Blowflies Model", xlab = "Days", ylab = "Adult Blowfly Counts")
points(x = blowfly.day, y = blowfly.data)
dev.off()

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


## Generating Data
blowfly.data <- matrix(blowfly.data, length(blowfly.data),1)
colnames(blowfly.data) <- c("y")
blowfly.data.d <- blowfly.data[times0 >= 20, , drop=FALSE]

# To get an initial estimate of the states we smooth the observed I component
# and set the other coefficients to zero.

# smooth the log I values
fdnames=list(NULL,c('y'),NULL)
DEfd0 <- smooth.basis(times0 ,(blowfly.data),fdPar(bbasis0,1,0.1))
DEfd.d <- smooth.basis(times.d, blowfly.data.d,fdPar(bbasis.d,1,0.1))
coefs0 <-  DEfd0$fd$coefs
coefs.d <- DEfd.d$fd$coefs
colnames(coefs0) <- colnames(coefs.d) <- c("y")
# set up the functional data object for the three variables
# plot the smooth plus data
## plotfit.fd(blowfly.data[,2],times0,DEfd0$fd)
DEfd0 <- fd(coefs0,bbasis0, fdnames)
DEfd.d <- fd(coefs.d,bbasis.d, fdnames)



#  list object betamore is now passed into LS.setup, too, in order to make
#  it available as a functional parameter defined by its three coefficients
#  run LS.setup
initBeta <- rep(0, 10)
initBeta[5:7] <- 1/5
initBeta[6] <- 3/5
blowfly.pars <- c(150 / 8, 3 / 8 , 1000)
names(blowfly.pars) <- c("c", "a", "N0")

##args <- commandArgs(TRUE)
## lambda <- 10^(-as.numeric(args[1])) * 10
lambda <- 1000

initBeta
blowfly.pars

dde.fit <- Profile.LS.sparse(blowfliesfn, blowfly.data.d, times.d, pars = blowfly.pars, beta = initBeta, coefs = coefs.d, basisvals = bbasis.d, lambda = lambda, in.meth='nlminb', delay = delay, basisvals0 = bbasis0, coefs0 = coefs0, nbeta = length(initBeta), ndelay = 1, tau = list(seq(3,12, by = 1)), control.out = list(method = "nnls", maxIter = 20, lambda.sparse = 0, echo = TRUE))

DEfd.fit <- fd(dde.fit$res$coefs, bbasis.d)
pdf(file = "blowfly-sim-fit.pdf", width = 7, height = 5)
plotfit.fd(blowfly.data.d,times.d,DEfd.fit, main = "Nicholson's Blowflies Model", xlab = "Days", ylab = "Adult Blowfly Counts")
dev.off()


lambda <- 1000
initBeta <- init.beta.zoom(dde.fit$res$beta[4:9])
dde.fit1 <- Profile.LS.sparse(blowfliesfn, blowfly.data.d, times.d, pars = dde.fit$res$pars, beta = initBeta, coefs = dde.fit$res$coefs, basisvals = bbasis.d, lambda = lambda, in.meth='nlminb', delay = delay, basisvals0 = bbasis0, coefs0 = coefs0, nbeta = length(initBeta), ndelay = 1, tau = list(seq(6,11, by = 0.5)), control.out = list(method = "nnls", maxIter = 20, lambda.sparse = 0, echo = TRUE))
DEfd.fit <- fd(dde.fit1$res$coefs, bbasis.d)
plotfit.fd(blowfly.data.d,times.d,DEfd.fit)

dde.fit2 <- Profile.LS.sparse(blowfliesfn, blowfly.data.d, times.d, pars = dde.fit$res$pars, beta = dde.fit$res$beta, coefs = coefs.d, basisvals = bbasis.d, lambda = lambda, in.meth='nlminb', delay = delay, basisvals0 = bbasis0, coefs0 = coefs0, nbeta = length(initBeta), ndelay = 1, tau = list(seq(3,12, by = 1)), control.out = list(method = "nnls", maxIter = 20, lambda.sparse = 0, echo = TRUE))



save(dde.fit2, dde.fit, dde.fit1, lambda, file = paste("blowfly-fit",lambda,".RData", sep=""))
