library(reshape2)
library(ggplot2)
library(penalized)
library(CollocInfer)
library(limSolve)

source("./R/tv-delay.R")
source("./R/mDTVSIRfn.R")
source("./R/sparse.R")
source("./R/LS.sparse.R")



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
        if(i == 2 && j == 2) {
        }else{
            load(paste("tv-fit", 100^(i-1)/1000, 100^(j-1)/1000, ".RData", sep=""))
            pars[[i]][[j]] <- tv.fit$res$pars
            beta[[i]][[j]] <- tv.fit$res$beta
            kappa[[i]][[j]] <- tv.fit$res$kappa
            fconv[[i]][j] <- tv.fit$res$f
        }
    }
}

i <- 3
j <- 1
load(paste("tv-fit", 100^(i-1)/1000, 100^(j-1)/1000, ".RData", sep=""))

DEfd.fit <- fd(tv.fit$res$coefs[,2, drop = FALSE],bbasis.d)
pdf()
plotfit.fd(mData.d[,2],times.d,DEfd.fit)
dev.off()



##################################################
## Parse 25 measles nnls fitting results
## 4 time varying coefficients
## Simulation script: measles-scr02.R
## Wed Oct  1 14:05:15 CDT 2014
## Commit: a5dc7d8e4c6c888789639092779a19c1e00be85b
##################################################

pars <- list()
beta <- list()
kappa <- list()
fconv <- list()

for(i in 1:5){
    beta[[i]] <- kappa[[i]]  <- list()
    pars[[i]] <- 0
    fconv[[i]] <- NA
    for(j in 1:5){
        if(i == 5 & j==5){}
        else{
            load(paste("tv-fit02", 100^(i-1)/1000, 100^(j-1)/1000, ".RData", sep=""))
            pars[[i]][j] <- tv.fit$res$pars
            beta[[i]][[j]] <- tv.fit$res$beta
            kappa[[i]][[j]] <- tv.fit$res$kappa
            fconv[[i]][j] <- tv.fit$res$f
        }
    }
}

##################################################
## Keep!!!
##################################################

i <- 4
j <- 1
load(paste("tv-fit02", 100^(i-1)/1000, 100^(j-1)/1000, ".RData", sep=""))

DEfd.fit <- fd(tv.fit$res$coefs[,2, drop = FALSE],bbasis.d)
#pdf()
plotfit.fd(mData.d[,2],times.d,DEfd.fit)
#dev.off()



##################################################
## Parse 25 measles nnls fitting results
## 4 time varying coefficients
## Simulation script: measles-scr03.R
## Tue Nov 11 16:24:08 CST 2014
## Commit: 622a8e1561fb2295cc1e302105fa4f44fd67ebcc
##################################################
## Using data from 1948 to 1963
## 12 time varying coefficients, uniform initial values
## 7 possible lags of delay, 0.5 at first two.

pars <- list()
beta <- list()
kappa <- list()
fconv <- list()

for(i in 1:5){
    beta[[i]] <- kappa[[i]]  <- list()
    pars[[i]] <- list()
    fconv[[i]] <- NA
    for(j in 1:5){
            tryCatch(
            {
                load(paste("mfit03-", 10^(i-1)/1000, 10^(j-1)/1000, ".RData", sep=""))
                pars[[i]][[j]] <- tv.fit$res$pars
                beta[[i]][[j]] <- tv.fit$res$beta
                kappa[[i]][[j]] <- tv.fit$res$kappa
                fconv[[i]][j] <- tv.fit$res$f
            },
                error= function(cond){
                    message(cond)
                    message(paste("\n", "Error at i,j being", i, j))
                    return(NA)
                }
                )
    }
}

i <- 2
j <- 2
load(paste("mfit03-", 10^(i-1)/1000, 10^(j-1)/1000, ".RData", sep=""))


i <- 4
j <- 2
load(paste("mfit03-", 10^(i-1)/1000, 10^(j-1)/1000, ".RData", sep=""))
DEfd.fit <- fd(tv.fit$res$coefs[,2, drop = FALSE],bbasis.d)
pdf()
plotfit.fd(mData.d[,2],times.d,DEfd.fit)
dev.off()


##################################################
## Parse 25 measles nnls fitting results
## 4 time varying coefficients
## Simulation script: measles-scr04.R
## Sun Nov 16 22:18:01 CST 2014
## Commit: bc69de0407396cb92170a235d6fc966ffd499eff
##################################################
## Using data from 1948 to 1963
## 12 time varying coefficients, uniform initial values
## 7 possible lags of delay, 0.5 at first two.

pars <- list()
beta <- list()
kappa <- list()
fconv <- list()

for(i in 1:4){
    beta[[i]] <- kappa[[i]]  <- list()
    pars[[i]] <- list()
    fconv[[i]] <- NA
    for(j in 1:4){
            tryCatch(
            {
                load(paste("mfit04-", 10^(i-1)/1000, 10^(j-1)/1000, ".RData", sep=""))
                pars[[i]][[j]] <- tv.fit$res$pars
                beta[[i]][[j]] <- tv.fit$res$beta
                kappa[[i]][[j]] <- tv.fit$res$kappa
                fconv[[i]][j] <- tv.fit$res$f
            },
                error= function(cond){
                    message(cond)
                    message(paste("\n", "Error at i,j being", i, j))
                    return(NA)
                }
                )
        }
}

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

i <- 1
j <- 1
load(paste("mfit04-", 10^(i-1)/1000, 10^(j-1)/1000, ".RData", sep=""))
DEfd.fit <- fd(tv.fit$res$coefs[,2, drop = FALSE],bbasis.d)
pdf()
plotfit.fd(mData.d[,2],times.d,DEfd.fit)
dev.off()

