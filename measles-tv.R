source("./R/timevar.R")
source("./R/mTVSIRfn.R")

library(penalized)
library(CollocInfer)
library(deSolve)
library(nnls)

mDf <- read.csv("./Data/meas_ca_on__1939-89_wk.csv", skip = 3)
bDf <- read.csv("./Data/bth_ca_on__1921-2002_mn.csv", skip = 5)
plot(x = mDf$numdate[mDf$numdate > 1955 & mDf$numdate < 1960], y = mDf$cases[mDf$numdate > 1955 & mDf$numdate < 1960],type = "l")

mTimes <- mDf$numdate[mDf$numdate > 1958 & mDf$numdate < 1963]
mI <- mDf$cases[mDf$numdate > 1958 & mDf$numdate < 1963]
tmpMonth <- mDf$month[mDf$numdate > 1958 & mDf$numdate < 1963]
mB <- rep(0, length(tmpMonth))

for(i in 1:length(tmpMonth)){
    mB[i] <- bDf[which(bDf$year == floor(mTimes[i]) & bDf$month == tmpMonth[i]),3 ] * 12 * 0.384 * (671 + 772) / 11410
}

mTimes <- mTimes - 1958
rr     = c(0,round(max(mTimes)))       #  the range of observations times
knots  = seq(rr[1],rr[2],2/52)  #  knots at 52 equally spaced values
norder = 3                      #  the order of the B-spline basis functions,
                                #  in this case piece-wise quadratic
nbasis = length(knots)+norder-2 #  the number of basis functions
bbasis <- create.bspline.basis(range=rr, norder=norder, nbasis=nbasis,breaks=knots)

mData <- matrix(NA, length(mI),2)
mData[,2] <- mI
colnames(mData) <- c("S" , "I")

# To get an initial estimate of the states we smooth the observed I component
# and set the other coefficients to zero.
fdnames=list(NULL,c('S', 'I'),NULL)
DEfd <- smooth.basis(mTimes ,(mData[,2]),fdPar(bbasis, 1, 0.1))
coefs <- cbind(matrix(800000,bbasis$nbasis,1), DEfd$fd$coefs)
colnames(coefs) = c("S","I")

procTimes <- c(0, seq(1/52, 5 - 1/52, by = 2/52), 5)

for(i in 1:length(procTimes)){
    month <- round((procTimes[i] - floor(procTimes[i])) * 12)
    if(month == 0)
        month <- 1
    procB[i] <- bDf[which(bDf$year == (floor(procTimes[i])+1958) & bDf$month == month),3 ] * 0.384 * 12
}

## * 12 * 0.384 * (671 + 772) / 11410 comes from:
## http://en.wikipedia.org/wiki/Demographics_of_Ontario
## \times 12 gives yearly birth rate
## \times 0.384 because of the Percentage of National Population of Ontario province
## \times (671 + 772) / 11410 because of age group

mPars <- c(20, 0.025)
names(mPars) <- c("gamma", "alpha")
mKappa <- rep(1e-3, 12)
names(mKappa) <- c("k1", "k2", "k3","k4","k5","k6","k7","k8","k9","k10","k11", "k12")

## debug(Profile.LS.tv)
tv.fit <- Profile.LS.tv(mTVSIRfn, mData, mTimes, pars = mPars, kappa = mKappa, coefs = coefs, basisvals = bbasis, lambda = c(1,10000), more = list(b = procB), in.meth='nlminb', control.out = list(method = "nnls", maxIter = 10, lambda.sparse = 0, echo = TRUE))

DEfd.fit <- fd(tv.fit$res$coefs[,2, drop = FALSE], bbasis)
plotfit.fd(mData[,2],mTimes,DEfd.fit)
DEfd.fit <- fd(tv.fit$res$coefs[,1, drop = FALSE] ,bbasis)
plotfit.fd(y = mData[,2] , argvals = mTimes, fdobj = DEfd.fit)

tv.fit$res$pars
tv.fit$res$kappa
