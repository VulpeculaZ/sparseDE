library(gpDDE)
library(spam)
load("data-blowfly-500.RData")

cov.one <- function(l){
     pars.hat <- l$pars
     beta.hat <- l$beta
     blowfly.data <- l$blowfly.data
     blowfly.data.d <- blowfly.data[times0 >= 20]
     blowfly.data <- matrix(blowfly.data, length(blowfly.data),1)
     blowfly.data.d <- matrix(blowfly.data.d, length(blowfly.data.d),1)
     coefs <- l$coefs
     DEfd0 <- smooth.basis(times0, blowfly.data,fdPar(bbasis0,1,0.1))
     coefs0 <-  DEfd0$fd$coefs
     Covar <- ProfileSSE.covariance.DDE(pars = pars.hat, beta = beta.hat, active = NULL, fn = blowfliesfn, data = blowfly.data.d, times = times.d,  coefs = coefs, basisvals = bbasis.d, lambda = lambda, in.meth='nlminb', basisvals0 = bbasis0, coefs0 = coefs0, nbeta = length(beta.hat), ndelay = 1, tau = tau)
     allpars <- c(pars.hat, beta.hat)
     cover <- ((allpars + 1.96 * sqrt(diag(Covar))) >= c(pars.true, beta.true)) & ((allpars - 1.96 * sqrt(diag(Covar))) <= c(pars.true, beta.true))
     return(list(Covar = Covar, coverage = cover))
}

## colnames(beta.hat) <- paste("beta", 1:10, sep = ".")
## beta.df <- melt(as.data.frame(beta.hat))
## pdf(file = "blowfly-250.pdf", width = 9,height = 5)
## ggplot(beta.df ,aes(x = variable,y = value))  + geom_boxplot()
## dev.off()
## Confidence interval estimation

blowfly.day <- seq(0,175, 0.5)
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
bfdPar0 = fdPar(bbasis0,lambda=1,int2Lfd(1))
bfdPar.d <- fdPar(bbasis.d,lambda=1,int2Lfd(1))
fdnames=list(NULL,c('y'),NULL)
lambda <- 10000
tau <- list(seq(5.5,10,0.5))
beta.true <- rep(0, 10)
beta.true[6] <- 1
pars.true <- c(150 / 8, 8 / 8 , 1000)
source("./R/make.blowfly.R")
blowfliesfn <- make.blowfly()
args <- commandArgs(TRUE)
i <- as.numeric(args[1])

nnls.res.all <- list()
load(paste("blowfly-nnls-500-true", i, ".RData", sep=""))
for(j in 1:length(nnls.res)){
    nnls.res.all[[j]] <- c(nnls.res[[j]], blowfly.data = list(data.res[[i*25 + j]]$blowfly.data))
}

cov.all <- list()
for(j in 1:length(nnls.res.all)){
    print(j)
    try(cov.all[[j]] <- cov.one(nnls.res.all[[j]]))
}
save(cov.all, file = paste("parse-cov-500-true", i, ".RData", sep = ""))


## > colMeans(coverage)
##       c        a       N0  beta1.1  beta1.2  beta1.3  beta1.4  beta1.5
##   1.000    0.988    0.988    1.000    1.000    1.000    1.000    1.000
## beta1.6  beta1.7  beta1.8  beta1.9 beta1.10
##   1.000    1.000    1.000    1.000    1.000
