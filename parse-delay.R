##################################################
## Parse 100 lasso fitting results of two seperate delays
## Lasso are used in every outer iteration
## Simulation script:
## sim.2d.lasso.R
## Commit: df114f8982bc1c11c228a4de5ca507198df08e66
## Tue Apr  1 15:55:22 CDT 2014
##################################################

library(CollocInfer)

times <- seq(0, 25, by = 0.1)
times.d <- knots.d <- times[times >= 5]
norder = 3
nbasis.d = length(knots.d) + norder - 2
range.d <- range(knots.d)
basis.d <- create.bspline.basis(range=range.d, nbasis=nbasis.d, norder=norder, breaks=knots.d)
bvals <- eval.basis(times.d, basis.d, 0)
bic <- matrix(NA, 50, 100)

for(i in 1:50){
    filenamei <- paste("lasso2d", i,".RData", sep = "")
    load(filenamei)
    if(length(sim.res) != 100)
        cat("Iterition prematurally terminated:",i, "\n")
    for(j in 1:length(sim.res)){
        devals = as.matrix(bvals%*%sim.res[[j]]$ncoefs)
        f <- sim.res[[j]]$data - devals
        sd.pen <- sd(f)
        ll.pen <- - sum(f^2) / (sd.pen^2) / 2 - length(f) * log(sd.pen)
        bic[i,j] <- -2 * ll.pen + (sum(sim.res[[j]]$res$beta > 0) + length(sim.res[[j]]$res$pars)) * log(length(f))
    }
}

minbic <- apply(bic, 2, min)
posminbic <- apply(sweep(bic, 2 , minbic, "=="),2,which)


pars.true <- 0.5
beta.true <- rep(0,16)
beta.true[c(4,7)] <- 2
betaMat <- matrix(NA, 100, 16)
pars <- c()
fdp.beta <- fnp.beta <- 0
for(i in 1:100){
    filenamei <- paste("lasso2d", posminbic[i],".RData", sep = "")
    load(filenamei)
    betaMat[i,] <- sim.res[[i]]$res$beta
    pars <- c(pars, sim.res[[i]]$res$pars)
    fdp.beta <-  fdp.beta + sum(sim.res[[i]]$res$beta[-c(4,7)] != 0) / sum(beta.true == 0)
    fnp.beta <- fnp.beta + sum(sim.res[[i]]$res$beta[c(4,7)] == 0) / 2
}

fdp.beta / 100
## 0.1857143
fnp.beta / 100
## 0.635

sum(colMeans(betaMat))
mean(pars)


##################################################
## Parse 500 nnls fitting results of two adjecent delays, sd = 0.01
## Simulation script:
## nnls-16d-2dadj.R
## Commit: df114f8982bc1c11c228a4de5ca507198df08e66
## Tue Apr  1 11:50:07 CDT 2014
##################################################

pars.hat <- beta.hat <- c()
beta.true <- rep(0, 16)
beta.true[c(7,8)] <- 2
pars.true <- 0.5

for(i in 0:19){
    load(paste("nnls-2dadj-sd01-", i, ".RData", sep=""))
    for(j in 1:25){
        pars.hat <- c(pars.hat, nnls.res[[j]]$pars)
        beta.hat <- rbind(beta.hat, nnls.res[[j]]$beta)
    }
}

cat("nnls")
fdp <- sum(beta.hat[,-c(7,8)] != 0) / sum(beta.hat != 0)
print(fdp)
fnp <- sum(beta.hat[,c(7,8)] == 0) / sum(beta.hat == 0)
print(fnp)

library(reshape2)
library(ggplot2)
colnames(beta.hat) <- round(seq(0,5, length.out = 16), digits =2)
beta.df <- melt(as.data.frame(beta.hat))
pdf(file = "nnls-16d-adj.pdf", width = 9,height = 5)
ggplot(beta.df ,aes(x = variable,y = value))  + geom_boxplot() + ylab(expression(hat(beta)[GP])) + xlab(expression(tau))
dev.off()


## Covariance:
Covar.mat <- c()
for(i in 0:19){
    load(paste("nnls-2dadj-sd02-", i, ".RData", sep=""))
    for(j in 1:25){
        if(!is.na(sim.Covar[[j]][1])){
            Covar.mat <- rbind(Covar.mat, diag(sim.Covar[[j]]))
        }else{
            print(i * 20 + j)
        }
    }
}

## Coverage of 95% CI:
cover.par <- ((pars.hat + 1.96 * sqrt(Covar.mat[,1])) >= 0.5) & ((pars.hat - 1.96 * sqrt(Covar.mat[,1])) <= 0.5)
cover.beta7 <- ((beta.hat[,7] + 1.96 * sqrt(Covar.mat[,8])) >= 2) & ((beta.hat[,7] - 1.96 * sqrt(Covar.mat[,8])) <= 2)
cover.beta8 <- ((beta.hat[,8] + 1.96 * sqrt(Covar.mat[,9])) >= 2) & ((beta.hat[,8] - 1.96 * sqrt(Covar.mat[,9])) <= 2)
sum(cover.beta7) / 500
sum(cover.beta8) / 500


##################################################
## Parse 500 penalized fitting results of two adjecent delays, sd = 0.01
## Penalization is done at the last iteration
## Simulation script:
## pen-16d-2dadj.R
## Commit: 471e1186e954baed9300525f4de0a657bf994793
## Thu Apr  3 10:49:45 CDT 2014
##################################################

pars.hat <- beta.hat <- c()
beta.true <- rep(0, 16)
beta.true[c(7,8)] <- 2
pars.true <- 0.5

for(i in 0:19){
    load(paste("pen-2dadj-sd01-", i, ".RData", sep=""))
    for(j in 1:25){
        pars.hat <- c(pars.hat, pen.res[[j]]$pars)
        beta.hat <- rbind(beta.hat, pen.res[[j]]$beta)
    }
}


cat("pen \n")
fdp <- sum(beta.hat[,-c(7,8)] != 0) / sum(beta.hat != 0)
print(fdp)
fnp <- sum(beta.hat[,c(7,8)] == 0) / sum(beta.hat == 0)
print(fnp)

colMeans(beta.hat)
sum(colMeans(beta.hat))
mean(pars.hat)
sd(pars.hat)

#library(reshape2)
#library(ggplot2)
colnames(beta.hat) <- round(seq(0,5, length.out = 16), digits =2)
beta.df <- melt(as.data.frame(beta.hat))
pdf(file = "pen-16d-adj.pdf", width = 9,height = 5)
ggplot(beta.df ,aes(x = variable,y = value))  + geom_boxplot() + ylab(expression(hat(beta)[GP-LASSO])) + xlab(expression(tau))
dev.off()


##################################################
## Parse 500 adaptive lasso fitting results of two adjecent delays, sd = 0.01
## Simulation script:
## al-16d-2dadj.R
## Commit: f771c451e6e3bbf8cd243346dd9044f84c533e3d
## Wed Apr  2 13:58:10 CDT 2014
##################################################

pars.hat <- beta.hat <- c()
beta.true <- rep(0, 16)
beta.true[c(7,8)] <- 2
pars.true <- 0.5

for(i in 0:19){
    load(paste("al-2dadj-sd01-", i, ".RData", sep=""))
    for(j in 1:25){
        pars.hat <- c(pars.hat, al.res[[j]]$pars)
        beta.hat <- rbind(beta.hat, al.res[[j]]$beta)
    }
}

cat("Adaptive LASSO \n")
fdp <- sum(beta.hat[,-c(7,8)] != 0) / sum(beta.hat != 0)
print(fdp)
fnp <- sum(beta.hat[,c(7,8)] == 0) / sum(beta.hat == 0)
print(fnp)

colMeans(beta.hat)
sum(colMeans(beta.hat))
mean(pars.hat)
sd(pars.hat)

#library(reshape2)
#library(ggplot2)
colnames(beta.hat) <- round(seq(0,5, length.out = 16), digits =2)
beta.df <- melt(as.data.frame(beta.hat))
pdf(file = "al-16d-adj.pdf", width = 9,height = 5)
ggplot(beta.df ,aes(x = variable,y = value))  + geom_boxplot() + ylab(expression(hat(beta)[GP-Adaptive-LASSO])) + xlab(expression(tau))
dev.off()

##################################################
## Parse 500 nnls fitting results of two adjecent delays, sd = 0.02
## Simulation script:
## nnls-16d-2dadj.R
## Commit: 9969ab5
## Sun Oct 25 11:41:52 CDT 2015
##################################################

pars.hat <- beta.hat <- c()
beta.true <- rep(0, 16)
beta.true[c(7,8)] <- 2
pars.true <- 0.5

for(i in 0:19){
    load(paste("nnls-2dadj-sd02-", i, ".RData", sep=""))
    for(j in 1:25){
        pars.hat <- c(pars.hat, nnls.res[[j]]$pars)
        beta.hat <- rbind(beta.hat, nnls.res[[j]]$beta)
    }
}

cat("nnls")
fdp <- sum(beta.hat[,-c(7,8)] != 0) / sum(beta.hat != 0)
print(fdp)
fnp <- sum(beta.hat[,c(7,8)] == 0) / sum(beta.hat == 0)
print(fnp)

colMeans(beta.hat)
sum(colMeans(beta.hat))
mean(pars.hat)
sd(pars.hat)
sd(beta.hat[,7])
sd(beta.hat[,8])
library(reshape2)
library(ggplot2)
colnames(beta.hat) <- round(seq(0,5, length.out = 16), digits =2)
beta.df <- melt(as.data.frame(beta.hat))
pdf(file = "nnls-16d-adj-sd02.pdf", width = 9,height = 5)
ggplot(beta.df ,aes(x = variable,y = value))  + geom_boxplot() + ylab(expression(hat(beta)[GP])) + xlab(expression(tau))
dev.off()



## Covariance:
Covar.mat <- c()
for(i in 0:19){
    load(paste("nnls-2dadj-sd02-", i, ".RData", sep=""))
    for(j in 1:25){
        if(!is.na(sim.Covar[[j]][1])){
            Covar.mat <- rbind(Covar.mat, diag(sim.Covar[[j]]))
        }else{
            print(i * 20 + j)
        }
    }
}

## Coverage of 95% CI:
pars.hat <- pars.hat[-c(68, 88)]
beta.hat <- beta.hat[-c(68, 88),]
cover.par <- ((pars.hat + 1.96 * sqrt(Covar.mat[,1])) >= 0.5) & ((pars.hat - 1.96 * sqrt(Covar.mat[,1])) <= 0.5)
cover.beta7 <- ((beta.hat[,7] + 1.96 * sqrt(Covar.mat[,8])) >= 2) & ((beta.hat[,7] - 1.96 * sqrt(Covar.mat[,8])) <= 2)
cover.beta8 <- ((beta.hat[,8] + 1.96 * sqrt(Covar.mat[,9])) >= 2) & ((beta.hat[,8] - 1.96 * sqrt(Covar.mat[,9])) <= 2)
sum(cover.beta7) / 498
sum(cover.beta8) / 498

##################################################
## Parse 500 lasso fitting results of two adjecent delays, sd = 0.02
## Simulation script:
## pen-16d-2dadj.R
## Commit: 3b9fcac
## Tue Oct 27 11:09:51 CDT 2015
##################################################

pars.hat <- beta.hat <- c()
beta.true <- rep(0, 16)
beta.true[c(7,8)] <- 2
pars.true <- 0.5

## Positive Lars:
for(i in 0:19){
    load(paste("pen-2dadj-sd02-", i, ".RData", sep=""))
    for(j in 1:25){
        pars.hat <- c(pars.hat, sim.res.lars[[j]]$pars)
        beta.hat <- rbind(beta.hat, sim.res.lars[[j]]$beta)
    }
}

cat("pen\n")
fdp <- sum(beta.hat[,-c(7,8)] != 0) / sum(beta.hat != 0)
print(fdp)
fnp <- sum(beta.hat[,c(7,8)] == 0) / sum(beta.hat == 0)
print(fnp)

colMeans(beta.hat)
sum(colMeans(beta.hat))
mean(pars.hat)
sd(pars.hat)

#library(reshape2)
#library(ggplot2)

colnames(beta.hat) <- round(seq(0,5, length.out = 16), digits =2)
beta.df <- melt(as.data.frame(beta.hat))
pdf(file = "pen-16d-adj-sd02.pdf", width = 9,height = 5)
ggplot(beta.df ,aes(x = variable,y = value))  + geom_boxplot() + ylab(expression(hat(beta)[GP-LASSO])) + xlab(expression(tau))
dev.off()


##################################################
## Parse 500 adaptive-lasso fitting results of two adjecent delays, sd = 0.02
## Simulation script:
## pen-16d-2dadj.R
## Commit: 3b9fcac
## Tue Oct 27 13:37:51 CDT 2015
##################################################

pars.hat <- beta.hat <- c()
beta.true <- rep(0, 16)
beta.true[c(7,8)] <- 2
pars.true <- 0.5

for(i in 0:19){
    load(paste("pen-2dadj-sd02-", i, ".RData", sep=""))
    for(j in 1:25){
        pars.hat <- c(pars.hat, sim.res.adlars[[j]]$pars)
        beta.hat <- rbind(beta.hat, sim.res.adlars[[j]]$beta)
    }
}

cat("Adaptive LASSO\n")
fdp <- sum(beta.hat[,-c(7,8)] != 0) / sum(beta.hat != 0)
print(fdp)
fnp <- sum(beta.hat[,c(7,8)] == 0) / length(beta.hat == 0)
print(fnp)

colMeans(beta.hat)
sum(colMeans(beta.hat))
mean(pars.hat)
sd(pars.hat)

#library(reshape2)
#library(ggplot2)
colnames(beta.hat) <- round(seq(0,5, length.out = 16), digits =2)
beta.df <- melt(as.data.frame(beta.hat))
pdf(file = "al-16d-adj-sd02.pdf", width = 9,height = 5)
ggplot(beta.df ,aes(x = variable,y = value))  + geom_boxplot() + ylab(expression(hat(beta)[GP-Adaptive-LASSO])) + xlab(expression(tau))
dev.off()



##################################################
## Parse 500 nnls fitting results of blowfly simulation, sd = 500
## Simulation script:
## blowfly-sim-batch.R
## Commit: fbde97e6f63178e3fcabf5d8e50de64b1a2a2b0b
## Tue Mar  3 21:24:36 CST 2015
##################################################

pars.hat <- beta.hat <- c()
beta.true <- rep(0, 10)
beta.true[6] <- 1
pars.true <- c(150 / 8, 8 / 8 , 1000)
Covar <- list()
coverage <- c()

for(i in 0:19){
    load(paste("blowfly-nnls-500-", i, ".RData", sep=""))
    for(j in 1:25){
        pars.hat <- rbind(pars.hat, nnls.res[[j]]$pars)
        beta.hat <- rbind(beta.hat, nnls.res[[j]]$beta)
    }
}

cat("nnls blowfly\n")
fdp <- sum(beta.hat[,-6] != 0) / length(beta.hat[,-6])
print(fdp)
fnp <- sum(beta.hat[,6] == 0) / length(beta.hat[,6])
print(fnp)
colMeans(beta.hat)
colMeans(pars.hat)

library(reshape2)
library(ggplot2)

colnames(beta.hat) <- seq(5.5,10,0.5)
beta.df <- melt(as.data.frame(beta.hat))
pdf(file = "blowfly-500-gp.pdf", width = 9,height = 5)
ggplot(beta.df ,aes(x = variable,y = value)) + ylab(expression(hat(beta))) + xlab(expression(tau)) + geom_boxplot()
dev.off()
rm(list = ls())



## Confidence interval estimation
load("data-blowfly-500.RData")
library(CollocInfer)
library(MASS)
source("./R/blowflies.R")
source("./R/sparse.R")
source("./R/LS.sparse.R")
source("./R/tv-delay-cov.R")

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
lambda <- 1000
tau <- list(seq(5.5,10,0.5))
beta.true <- rep(0, 10)
beta.true[6] <- 1
pars.true <- c(150 / 8, 8 / 8 , 1000)

Covar <- list()
coverage <- c()
blowfliesfn <- make.blowfly()

load("blowfly-nnls-500-0.RData")
DEfd.fit <- fd(nnls.res[[1]]$coefs, bbasis.d)
blowfly.data <- data.res[[1]]$blowfly.data
blowfly.data.d <- blowfly.data[times0 >= 20]

pdf(file = "blowfly-sim-fit.pdf", width = 8, height = 5)
plotfit.fd(blowfly.data.d,times.d, DEfd.fit, xlab = "t (days)", ylab = "y", titles = NA, sub = NA)
dev.off()

## Do not run
## Takes a long time.
## The functionality is wrapped in cov.one() and can be parallelized.
for(i in 0:19){
    load(paste("blowfly-nnls-500-", i, ".RData", sep=""))
    for(j in 1:25){
        pars.hat <- nnls.res[[j]]$pars
        beta.hat <- nnls.res[[j]]$beta
        blowfly.data <- data.res[[i*25 + j]]$blowfly.data
        blowfly.data.d <- blowfly.data[times0 >= 20]
        blowfly.data <- matrix(blowfly.data, length(blowfly.data),1)
        blowfly.data.d <- matrix(blowfly.data.d, length(blowfly.data.d),1)
        coefs <- nnls.res[[j]]$coefs
        DEfd0 <- smooth.basis(times0, blowfly.data,fdPar(bbasis0,1,0.1))
        coefs0 <-  DEfd0$fd$coefs
        Covar[[i*25 + j]] <- ProfileSSE.covariance.delay(pars = pars.hat, beta = beta.hat, active = NULL, fn = blowfliesfn, data = blowfly.data.d, times = times.d,  coefs = coefs, basisvals = bbasis.d, lambda = lambda, in.meth='nlminb', delay = delay, basisvals0 = bbasis0, coefs0 = coefs0, nbeta = length(beta.hat), ndelay = 1, tau = tau)
        allpars <- c(pars.hat, beta.hat)
        cover <- ((allpars + 1.96 * sqrt(diag(Covar[[i*25 + j]]))) >= c(pars.true, beta.true)) & ((allpars - 1.96 * sqrt(diag(Covar[[i*25 + j]]))) <= c(pars.true, beta.true))
        coverage <- rbind(coverage, cover)
    }
}

save.image()
 ##       c        a       N0  beta1.1  beta1.2  beta1.3  beta1.4  beta1.5
 ##     452      409      389      499      499      500      500      499
 ## beta1.6  beta1.7  beta1.8  beta1.9 beta1.10
 ##     496      497      500      500      500

nnls.res.all <- list()
for(i in 0:19){
    load(paste("blowfly-nnls-500-", i, ".RData", sep=""))
    for(j in 1:25){
        nnls.res[[j]] <- c(nnls.res[[j]], blowfly.data = list(data.res[[i*25 + j]]$blowfly.data))
    }
    nnls.res.all <- c(nnls.res.all, nnls.res)
}

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
tmp <- cov.one(nnls.res.all[[1]])

library(parallel)
system.time(mclapply(nnls.res.all[1:50], cov.one, mc.cores = 25))
system.time(lapply(nnls.res.all[1:2], cov.one))
cov.all <- mclapply(nnls.res.all[1:50], cov.one, mc.cores = 25)

##################################################
## LASSO and ad-lars
##################################################

pars.hat <- beta.hat <- c()
beta.true <- rep(0, 10)
beta.true[6] <- 1
pars.true <- c(150 / 8, 8 / 8 , 1000)
Covar <- list()
coverage <- c()

for(i in 0:19){
    load(paste("blowfly-lasso-500-", i, ".RData", sep=""))
    for(j in 1:25){
        pars.hat <- rbind(pars.hat, sim.res.lars[[j]]$pars)
        beta.hat <- rbind(beta.hat, sim.res.lars[[j]]$beta)
    }
}

cat("nnls blowfly\n")
fdp <- sum(beta.hat[,-6] != 0) / length(beta.hat[,-6])
print(fdp)
fnp <- sum(beta.hat[,6] == 0) / length(beta.hat[,6])
print(fnp)
colMeans(beta.hat)
colMeans(pars.hat)

library(reshape2)
library(ggplot2)
colnames(beta.hat) <- seq(5.5,10,0.5)
beta.df <- melt(as.data.frame(beta.hat))
pdf(file = "blowfly-500-lasso.pdf", width = 9,height = 5)
ggplot(beta.df ,aes(x = variable,y = value)) + ylab(expression(hat(beta))) + xlab(expression(tau)) + geom_boxplot()
dev.off()
rm(list = ls())

##################################################
## Parse 500 nnls fitting results of blowfly simulation, sd = 1000
## Simulation script:
## blowfly-sim-batch.R
## Commit: 80de87543838ad14ab57b8db37d7a23f9756adf3
## Sun Mar 15 20:01:49 CDT 2015
##################################################

pars.hat <- beta.hat <- c()
beta.true <- rep(0, 10)
beta.true[6] <- 1
pars.true <- c(150 / 8, 8 / 8 , 1000)
Covar <- list()
coverage <- c()

for(i in 0:19){
    load(paste("blowfly-nnls-1000-", i, ".RData", sep=""))
    print(length(nnls.res))
    for(j in 1:length(nnls.res)){
        pars.hat <- rbind(pars.hat, nnls.res[[j]]$pars)
        beta.hat <- rbind(beta.hat, nnls.res[[j]]$beta)
    }
}

cat("nnls blowfly\n")
fdp <- sum(beta.hat[,-6] != 0) / length(beta.hat[,-6])
print(fdp)
fnp <- sum(beta.hat[,6] == 0) / length(beta.hat[,-6])
print(fnp)
colMeans(beta.hat)
colMeans(pars.hat)

library(reshape2)
library(ggplot2)
colnames(beta.hat) <- paste("beta", 1:10, sep = ".")
beta.df <- melt(as.data.frame(beta.hat))
ggplot(beta.df ,aes(x = variable,y = value))  + geom_boxplot()

## Confidence interval estimation
load("data-blowfly-1000.RData")
library(CollocInfer)
library(MASS)
source("./R/blowflies.R")
source("./R/sparse.R")
source("./R/LS.sparse.R")
source("./R/tv-delay-cov.R")

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
lambda <- 1000
tau <- list(seq(5.5,10,0.5))
beta.true <- rep(0, 10)
beta.true[6] <- 1
pars.true <- c(150 / 8, 8 / 8 , 1000)


nnls.res.all <- list()
for(i in 0:19){
    load(paste("blowfly-nnls-1000-", i, ".RData", sep=""))
    for(j in 1:length(nnls.res)){
        nnls.res[[j]] <- c(nnls.res[[j]], blowfly.data = list(data.res[[i*25 + j]]$blowfly.data))
    }
    nnls.res.all <- c(nnls.res.all, nnls.res)
}

library(parallel)
system.time(cov.all <- mclapply(nnls.res.all, cov.one, mc.cores = 30, mc.preschedule = FALSE))

load("parse-cov-1000.RData")
coverage <- c()
for(i in 1:length(cov.all)){
    try(coverage <- rbind(coverage, cov.all[[i]]$coverage))
}
colMeans(coverage)

##################################################
## Parse 500 nnls fitting results of blowfly simulation, sd = 250
## Simulation script:
## blowfly-sim-batch.R
## Commit: 5c0aa56905a66c6b6ff585fddeef0a42324dd43b
## Tue Mar 17 15:06:19 CDT 2015
##################################################

pars.hat <- beta.hat <- c()
beta.true <- rep(0, 10)
beta.true[6] <- 1
pars.true <- c(150 / 8, 8 / 8 , 1000)
Covar <- list()
coverage <- c()

for(i in 0:19){
    load(paste("blowfly-nnls-250-", i, ".RData", sep=""))
    print(length(nnls.res))
    for(j in 1:length(nnls.res)){
        pars.hat <- rbind(pars.hat, nnls.res[[j]]$pars)
        beta.hat <- rbind(beta.hat, nnls.res[[j]]$beta)
    }
}

cat("nnls blowfly\n")
fdp <- sum(beta.hat[,-6] != 0) / length(beta.hat[,-6])
print(fdp)
fnp <- sum(beta.hat[,6] == 0) / length(beta.hat[,-6])
print(fnp)
colMeans(beta.hat)
colMeans(pars.hat)

for(i in 1:3){
print(sd(pars.hat[,i]))
}

library(reshape2)
library(ggplot2)
colnames(beta.hat) <- seq(5.5,10,0.5)
beta.df <- melt(as.data.frame(beta.hat))
pdf(file = "blowfly-250.pdf", width = 9,height = 5)
ggplot(beta.df ,aes(x = variable,y = value)) + ylab(expression(hat(beta))) + xlab(expression(tau)) + geom_boxplot()
dev.off()

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
lambda <- 1000
tau <- list(seq(5.5,10,0.5))

load("data-blowfly-250.RData")
nnls.res.all <- list()
for(i in 0:19){
    load(paste("blowfly-nnls-250-", i, ".RData", sep=""))
    for(j in 1:length(nnls.res)){
        nnls.res[[j]] <- c(nnls.res[[j]], blowfly.data = list(data.res[[i*25 + j]]$blowfly.data))
    }
    nnls.res.all <- c(nnls.res.all, nnls.res)
}

tmp <- cov.one(nnls.res.all[[1]])

library(parallel)
cov.all <- mclapply(nnls.res.all[1:2], cov.one, mc.preschedule = FALSE)
save.image()

load("parse-cov-250.RData")
coverage <- c()
for(i in 1:length(cov.all)){
    try(coverage <- rbind(coverage, cov.all[[i]]$coverage))
}


## > colMeans(coverage)
##       c        a       N0  beta1.1  beta1.2  beta1.3  beta1.4  beta1.5
##   1.000    0.988    0.988    1.000    1.000    1.000    1.000    1.000
## beta1.6  beta1.7  beta1.8  beta1.9 beta1.10
##   1.000    1.000    1.000    1.000    1.000

fdpfnp <- function(fd, fn){
    fdp <- fd * 9 / (fd * 9 + 1 - fn)
    fnp <- fn / ((1 - fd) * 9 + fn)
    c(fdp, fnp)
}



##################################################
## Parse 500 nnls fitting results of blowfly simulation, sd = 250
## Simulation script:
## blowfly-sim-batch.R
## Commit: 5581fc49055bfb17fdbfc722d41b0294c5d17464
##
##################################################

pars.hat <- beta.hat <- c()
beta.true <- rep(0, 10)
beta.true[6] <- 1
pars.true <- c(150 / 8, 8 / 8 , 1000)
Covar <- list()
coverage <- c()

for(i in 0:19){
    load(paste("blowfly-nnls-250-dist", i, ".RData", sep=""))
    print(length(nnls.res))
    for(j in 1:length(nnls.res)){
        pars.hat <- rbind(pars.hat, nnls.res[[j]]$pars)
        beta.hat <- rbind(beta.hat, nnls.res[[j]]$beta)
    }
}

cat("nnls blowfly\n")
fdp <- sum(beta.hat[,-6] != 0) / length(beta.hat[,-6])
print(fdp)
fnp <- sum(beta.hat[,6] == 0) / length(beta.hat[,-6])
print(fnp)
colMeans(beta.hat)
colMeans(pars.hat)

library(reshape2)
library(ggplot2)
colnames(beta.hat) <- paste("beta", 1:10, sep = ".")
beta.df <- melt(as.data.frame(beta.hat))
pdf(file = "blowfly-250.pdf", width = 9,height = 5)
ggplot(beta.df ,aes(x = variable,y = value))  + geom_boxplot()
dev.off()
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
lambda <- 1000
tau <- list(seq(5.5,10,0.5))

load("data-blowfly-250.RData")
nnls.res.all <- list()
for(i in 0:19){
    load(paste("blowfly-nnls-250-", i, ".RData", sep=""))
    for(j in 1:length(nnls.res)){
        nnls.res[[j]] <- c(nnls.res[[j]], blowfly.data = list(data.res[[i*25 + j]]$blowfly.data))
    }
    nnls.res.all <- c(nnls.res.all, nnls.res)
}

tmp <- cov.one(nnls.res.all[[1]])

library(parallel)
cov.all <- mclapply(nnls.res.all[1:2], cov.one, mc.preschedule = FALSE)
save.image()

load("parse-cov-250.RData")
coverage <- c()
for(i in 1:length(cov.all)){
    try(coverage <- rbind(coverage, cov.all[[i]]$coverage))
}


## > colMeans(coverage)
##       c        a       N0  beta1.1  beta1.2  beta1.3  beta1.4  beta1.5
##   1.000    0.988    0.988    1.000    1.000    1.000    1.000    1.000
## beta1.6  beta1.7  beta1.8  beta1.9 beta1.10
##   1.000    1.000    1.000    1.000    1.000

fdpfnp <- function(fd, fn){
    fdp <- fd * 9 / (fd * 9 + 1 - fn)
    fnp <- fn / ((1 - fd) * 9 + fn)
    c(fdp, fnp)
}



##################################################
## Parse 500 nnls fitting results of blowfly simulation, sd = 500
## Simulation script:
## blowfly-sim-batch.R
## Commit: fbde97e6f63178e3fcabf5d8e50de64b1a2a2b0b
##
##################################################

pars.hat <- beta.hat <- c()
beta.true <- rep(0, 10)
beta.true[6] <- 1
pars.true <- c(150 / 8, 8 / 8 , 1000)
Covar <- list()
coverage <- c()

for(i in 0:19){
    load(paste("blowfly-nnls-500-true", i, ".RData", sep=""))
    load(paste("parse-cov-500-true",i,".RData", sep = ""))
    for(j in 1:25){
        pars.hat <- rbind(pars.hat, nnls.res[[j]]$pars)
        beta.hat <- rbind(beta.hat, nnls.res[[j]]$beta)
        try(coverage <- rbind(coverage, cov.all[[j]]$coverage))
    }
}

cat("nnls blowfly\n")
fdp <- sum(beta.hat[,-6] != 0) / sum(beta.hat[,-6] != 0)
print(fdp)
fnp <- sum(beta.hat[,6] == 0) / sum(beta.hat[,6] == 0)
print(fnp)
colMeans(beta.hat)
colMeans(pars.hat)

library(reshape2)
library(ggplot2)

colnames(beta.hat) <- seq(5.5,10,0.5)
beta.df <- melt(as.data.frame(beta.hat))
pdf(file = "blowfly-500-gp-true.pdf", width = 9,height = 5)
ggplot(beta.df ,aes(x = variable,y = value)) + ylab(expression(hat(beta))) + xlab(expression(tau)) + geom_boxplot()
dev.off()
rm(list = ls())

colSums(coverage) / 500

load("parse-cov-250.RData")
coverage <- c()
for(i in 1:length(cov.all)){
    try(coverage <- rbind(coverage, cov.all[[i]]$coverage))
}

       c        a       N0  beta1.1  beta1.2  beta1.3  beta1.4  beta1.5
   0.954    0.914    0.900    0.994    1.000    1.000    0.996    1.000
 beta1.6  beta1.7  beta1.8  beta1.9 beta1.10
   0.980    0.992    0.994    1.000    0.998




##################################################
## Parse 500 nnls fitting results of blowfly simulation, sd = 250
## Simulated distributed delay
## Simulation script:
## blowfly-sim-batch.R
## Commit: e3aeaa60a380e119969c555e802583a6fa07556b
## Mon Jan 25 12:15:16 CST 2016
##################################################

pars.hat <- beta.hat <- c()
beta.true <- rep(0, 10)
beta.true[6] <- 1
pars.true <- c(150 / 8, 8 / 8 , 1000)
Covar <- list()
coverage <- c()

for(i in 0:19){
    load(paste("blowfly-nnls-250-dist", i, ".RData", sep=""))
    for(j in 1:25){
        pars.hat <- rbind(pars.hat, nnls.res[[j]]$pars)
        beta.hat <- rbind(beta.hat, nnls.res[[j]]$beta)
    }
}

pars.hat <- pars.hat[-77,]
cat("nnls blowfly\n")
fdp <- sum(beta.hat[,-6] != 0) / sum(beta.hat[,-6] != 0)
print(fdp)
fnp <- sum(beta.hat[,6] == 0) / sum(beta.hat[,6] == 0)
print(fnp)
colMeans(beta.hat)
colMeans(pars.hat)

for(i in 1:3){
print(sd(pars.hat[,i]))
}

library(reshape2)
library(ggplot2)

colnames(beta.hat) <- seq(5.5,10,0.5)
beta.df <- melt(as.data.frame(beta.hat))
pdf(file = "blowfly-250-gp-dist.pdf", width = 9,height = 5)
ggplot(beta.df ,aes(x = variable,y = value)) + ylab(expression(hat(beta))) + xlab(expression(tau)) + geom_boxplot()
dev.off()
rm(list = ls())


##################################################
## Parse 500 nnls fitting results of blowfly simulation, sd = 500
## Simulation script:
## blowfly-sim-batch.R
## Commit: 6a77de5ebae55ff4d6f710c7922514965a8970ce
## Mon Jan 25 13:03:25 CST 2016
##################################################

pars.hat <- beta.hat <- c()
beta.true <- rep(0, 10)
beta.true[6] <- 1
pars.true <- c(150 / 8, 8 / 8 , 1000)
Covar <- list()
coverage <- c()

for(i in 0:19){
    load(paste("blowfly-nnls-500-true", i, ".RData", sep=""))
    load(paste("parse-cov-500-true",i,".RData", sep = ""))
    for(j in 1:25){
        pars.hat <- rbind(pars.hat, nnls.res[[j]]$pars)
        beta.hat <- rbind(beta.hat, nnls.res[[j]]$beta)
        try(coverage <- rbind(coverage, cov.all[[j]]$coverage))
    }
}

cat("nnls blowfly\n")
fdp <- sum(beta.hat[,-6] != 0) / sum(beta.hat[,-6] != 0)
print(fdp)
fnp <- sum(beta.hat[,6] == 0) / sum(beta.hat[,6] == 0)
print(fnp)
colMeans(beta.hat)
colMeans(pars.hat)

library(reshape2)
library(ggplot2)

colnames(beta.hat) <- seq(5.5,10,0.5)
beta.df <- melt(as.data.frame(beta.hat))
pdf(file = "blowfly-500-gp-true.pdf", width = 9,height = 5)
ggplot(beta.df ,aes(x = variable,y = value)) + ylab(expression(hat(beta))) + xlab(expression(tau)) + geom_boxplot()
dev.off()
rm(list = ls())

colSums(coverage) / 500

load("parse-cov-250.RData")
coverage <- c()
for(i in 1:length(cov.all)){
    try(coverage <- rbind(coverage, cov.all[[i]]$coverage))
}



##################################################
## Parse 500 nnls fitting results of two adjecent delays, sd = 0.01
## Simulation script:
## nnls-16d-2dadj.R
## Commit: df114f8982bc1c11c228a4de5ca507198df08e66
## Tue Apr  1 11:50:07 CDT 2014
##################################################

pars.hat <- beta.hat <- c()
beta.true <- rep(2,2)
pars.true <- 0.5
for(i in 0:19){
    load(paste("nnls-2dadj-sd02-true-", i, ".RData", sep=""))
    for(j in 1:25){
        pars.hat <- c(pars.hat, nnls.res[[j]]$pars)
        beta.hat <- rbind(beta.hat, nnls.res[[j]]$beta)
    }
}


library(reshape2)
library(ggplot2)
colnames(beta.hat) <- round(seq(0,5, length.out = 16), digits =2)
beta.df <- melt(as.data.frame(beta.hat))
pdf(file = "nnls-16d-adj.pdf", width = 9,height = 5)
ggplot(beta.df ,aes(x = variable,y = value))  + geom_boxplot() + ylab(expression(hat(beta)[GP])) + xlab(expression(tau))
dev.off()


## Covariance:
Covar.mat <- c()
for(i in 0:19){
    load(paste("nnls-2dadj-sd02-true-", i, ".RData", sep=""))
    for(j in 1:25){
        if(!is.na(sim.Covar[[j]][1])){
            Covar.mat <- rbind(Covar.mat, diag(sim.Covar[[j]]))
        }else{
            print(i * 20 + j)
        }
    }
}

## Coverage of 95% CI:
cover.par <- ((pars.hat + 1.96 * sqrt(Covar.mat[,1])) >= 0.5) & ((pars.hat - 1.96 * sqrt(Covar.mat[,1])) <= 0.5)
cover.beta1 <- ((beta.hat[,1] + 1.96 * sqrt(Covar.mat[,2])) >= 2) & ((beta.hat[,1] - 1.96 * sqrt(Covar.mat[,2])) <= 2)
cover.beta2 <- ((beta.hat[,2] + 1.96 * sqrt(Covar.mat[,3])) >= 2) & ((beta.hat[,2] - 1.96 * sqrt(Covar.mat[,3])) <= 2)
sum(cover.par) / 500
sum(cover.beta1) / 500
sum(cover.beta2) / 500

sum(cover.beta2 == 0)




##################################################
## Parse 500 nnls fitting results of one adjecent delays, sd = 0.02
## Simulation script:
## nnls-16d-2dadj.R
## Commit: 00b66d98aa8d339573b1c22cd46b838e2d8e505e
## Sun Jan 31 17:29:16 CST 2016
##################################################

pars.hat <- beta.hat <- c()
beta.true <- c(2)
pars.true <- 0.5
for(i in 0:19){
    load(paste("nnls-2dadj-sd02-1delay-", i, ".RData", sep=""))
    for(j in 1:25){
        pars.hat <- c(pars.hat, nnls.res[[j]]$pars)
        beta.hat <- rbind(beta.hat, nnls.res[[j]]$beta)
    }
}


library(reshape2)
library(ggplot2)
colnames(beta.hat) <- round(seq(0,5, length.out = 16), digits =2)
beta.df <- melt(as.data.frame(beta.hat))
pdf(file = "nnls-16d-adj.pdf", width = 9,height = 5)
ggplot(beta.df ,aes(x = variable,y = value))  + geom_boxplot() + ylab(expression(hat(beta)[GP])) + xlab(expression(tau))
dev.off()


## Covariance:
Covar.mat <- c()
for(i in 0:19){
    load(paste("nnls-2dadj-sd02-1delay-", i, ".RData", sep=""))
    for(j in 1:25){
        if(!is.na(sim.Covar[[j]][1])){
            Covar.mat <- rbind(Covar.mat, diag(sim.Covar[[j]]))
        }else{
            print(i * 20 + j)
        }
    }
}

## Coverage of 95% CI:
cover.par <- ((pars.hat + 1.96 * sqrt(Covar.mat[,1])) >= 0.5) & ((pars.hat - 1.96 * sqrt(Covar.mat[,1])) <= 0.5)
cover.beta1 <- ((beta.hat[,1] + 1.96 * sqrt(Covar.mat[,2])) >= 2) & ((beta.hat[,1] - 1.96 * sqrt(Covar.mat[,2])) <= 2)
sum(cover.par) / 500
sum(cover.beta1) / 500

sum(cover.beta2 == 0)
