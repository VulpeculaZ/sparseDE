library(reshape2)
library(ggplot2)

source("./R/tv-delay.R")
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


##################################################
## Parse 25 measles nnls fitting results
## Simulation script: measles-scr.R
## Wed Oct  1 14:05:15 CDT 2014
## Commit: 6682ea6af7092acc28ae5fe2f5f0abb1ce984736
##################################################

pars <- list()
beta <- list()
kappa <- list()
fconv <- list()

for(i in 1:5){
    pars[[i]] <- beta[[i]] <- kappa[[i]]  <- list()
    fconv[[i]] <- NA
    for(j in 1:5){
        load(paste("tv-fit", 10^(i-1), 10^(j-1), ".RData", sep=""))
        pars[[i]][[j]] <- tv.fit$res$pars
        beta[[i]][[j]] <- tv.fit$res$beta
        kappa[[i]][[j]] <- tv.fit$res$kappa
        fconv[[i]][j] <- tv.fit$res$f
    }
}

DEfd.fit <- fd(tv.fit$res$coefs[,2, drop = FALSE],bbasis.d)
pdf()
plotfit.fd(mData.d[,2],times.d,DEfd.fit)
dev.off()
