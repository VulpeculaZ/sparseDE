source("./R/tv-delay.R")
source("./R/mDTVtrendSIRfn.R")
source("./R/mDTVSIRfn.R")
source("./R/sparse.R")
source("./R/LS.sparse.R")

library(penalized)
library(CollocInfer)
library(limSolve)


## Real data:
mDf <- read.csv("./Data/meas_ca_on__1939-89_wk.csv", skip = 3)
bDf <- read.csv("./Data/bth_ca_on__1921-2002_mn.csv", skip = 5)
## plot(x = mDf$numdate[mDf$numdate > 1955 & mDf$numdate < 1960], y = mDf$cases[mDf$numdate > 1955 & mDf$numdate < 1960],type = "l")

mTimes <- mDf$numdate[mDf$numdate > 1948 & mDf$numdate <= 1963]
mI <- mDf$cases[mDf$numdate > 1948 & mDf$numdate <= 1963]
tmpMonth <- mDf$month[mDf$numdate > 1948 & mDf$numdate <= 1963]
mB <- rep(0, length(tmpMonth))

## pdf(file="measles-48to63-1.pdf", height = 6, width=10)
## plot( mTimes, mI, xlab = "Years", ylab = "Number of Infected", main = "Measles Cases in Ontario Province, Canada", type = "l")
## dev.off()

bTimes <- seq(1948, 1963 - 1/12, by = 1/12)
bB <-  bDf[which(bDf$year <= 1962 & bDf$year >= 1948),3]
## pdf(file="measles-48to63-birth.pdf", height = 6, width=10)
## plot( bTimes, bB, xlab = "Years", ylab = "Number of New Borns", main = "Monthly Birthrate  in Ontario Province, Canada", type = "l")
## dev.off()


for(i in 1:length(tmpMonth)){
    mB[i] <- bDf[which(bDf$year == floor(mTimes[i]) & bDf$month == tmpMonth[i]),3 ] * 12
}


mTimes <- mTimes - 1948
rr     = c(0,round(max(mTimes)))       #  the range of observations times
knots  = seq(rr[1],rr[2],2/52)  #  knots at 52 equally spaced values
norder = 3                      #  the order of the B-spline basis functions,
                                #  in this case piece-wise quadratic
nbasis = length(knots) + norder - 2 #  the number of basis functions

#  set up the basis object

bbasis0 <- create.bspline.basis(range=rr, norder=norder, nbasis=nbasis,breaks=knots)
times0  <- mTimes
times.d  <- mTimes[mTimes >= 1]
knots.d <-  seq(1,rr[2],2/52)
nbasis.d = length(knots.d) + norder - 2
bbasis.d <- create.bspline.basis(range=c(1,rr[2]), norder=norder, nbasis=nbasis.d, breaks=knots.d)


## Generating Data
mData <- matrix(NA, length(mI),2)
mData[,2] <- mI
colnames(mData) <- c("S" , "I")
mData.d <- mData[mTimes >= 1,]

# To get an initial estimate of the states we smooth the observed I component
# and set the other coefficients to zero.

# smooth the log I values
fdnames=list(NULL,c('S', 'I'),NULL)
DEfd0 <- smooth.basis(times0 ,(mData[,2]),fdPar(bbasis0,1,0.1))
DEfd.d <- smooth.basis(times.d, (mData[,2])[times0 >= 1],fdPar(bbasis.d,1,0.1))
coefs0 <- cbind(matrix(80000,bbasis0$nbasis,1), DEfd0$fd$coefs)
coefs.d <- cbind(matrix(80000,bbasis.d$nbasis,1), DEfd.d$fd$coefs)
colnames(coefs0) <- colnames(coefs.d) <- c("S", "I")
# set up the functional data object for the three variables
# plot the smooth plus data
## plotfit.fd(mData[,2],times0,DEfd0$fd)
DEfd0 <- fd(coefs0,bbasis0, fdnames)
DEfd.d <- fd(coefs.d,bbasis.d, fdnames)


procTimes <- c(1, seq(1 + 1/52, 15 - 1/52, by = 2/52), 15)
procB <- vector(,length(procTimes))
for(i in 1:length(procTimes)){
    month <- round((procTimes[i] - floor(procTimes[i])) * 12)
    if(month == 0)
        month <- 1
    procB[i] <- bDf[which(bDf$year == (floor(procTimes[i])+1948) & bDf$month == month),3 ] * 12
}

#  list object betamore is now passed into LS.setup, too, in order to make
#  it available as a functional parameter defined by its three coefficients
#  run LS.setup

args <- commandArgs(TRUE)
lambda1 <- 10^(as.numeric(args[1]) %/% 4) / 1000
lambda2 <- 10^(as.numeric(args[1]) %% 4) / 1000

initBeta <- rep(0, 7)
initBeta[1:2] <- 0.5
mPars <- c(mean(procB) / mean(mData.d[,2]), 0, 0)
names(mPars) <- c("gamma", "pho0", "pho1")
mKappa <- rep(2e-3, 12)
#mKappa[c(6,7,8)] <- 1e-3
names(mKappa) <- c("k1", "k2", "k3","k4","k5","k6","k7","k8","k9","k10","k11", "k12")


coefsS <- init.unob.LS.tv.delay(mDTVSIRfn, mData.d, times.d, pars = mPars, kappa = mKappa, coefs = coefs.d, beta = initBeta, basisvals = bbasis.d, lambda = c(lambda1,lambda2), more = list(b = procB), in.meth='nlminb', control.out = list(method = "nnls", maxIter = 10, lambda.sparse = 0, echo = TRUE), delay = delay, basisvals0 = bbasis0, coefs0 = coefs0, nbeta = length(initBeta), ndelay = 2, tau = list(seq(0,6/52, 1/52)), unob = 1)
coefs.d[, 1] <- coefsS$coefficients



## debug(Profile.LS.tv)
tv.fit <- Profile.LS.tv.delay(mDTVSIRtrfn, mData.d, times.d, pars = mPars, kappa = mKappa, coefs = coefs.d, beta = initBeta, basisvals = bbasis.d, lambda = c(lambda1,lambda2), more = list(b = procB), in.meth='nlminb', control.out = list(method = "nnls", maxIter = 10, lambda.sparse = 0, echo = TRUE), delay = delay, basisvals0 = bbasis0, coefs0 = coefs0, nbeta = length(initBeta), ndelay = 2, tau = list(seq(0,6/52, 1/52)))

save(tv.fit, lambda1, lambda2, file = paste("mfit04-",lambda1,lambda2,".RData", sep=""))
