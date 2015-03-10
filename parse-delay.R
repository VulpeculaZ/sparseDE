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
fdp <- sum(beta.hat[,-c(7,8)] != 0) / length(beta.hat[,-c(7,8)])
print(fdp)
fnp <- sum(beta.hat[,c(7,8)] == 0) / length(beta.hat[,c(7,8)])
print(fnp)

library(reshape2)
library(ggplot2)
colnames(beta.hat) <- paste("beta", 1:16, sep = ".")
beta.df <- melt(as.data.frame(beta.hat))
pdf(file = "nnls-16d-adj.pdf", width = 9,height = 5)
ggplot(beta.df ,aes(x = variable,y = value))  + geom_boxplot()
dev.off()

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
fdp <- sum(beta.hat[,-c(7,8)] != 0) / length(beta.hat[,-c(7,8)])
print(fdp)
fnp <- sum(beta.hat[,c(7,8)] == 0) / length(beta.hat[,c(7,8)])
print(fnp)
colMeans(beta.hat)
sum(colMeans(beta.hat))
mean(pars.hat)
sd(pars.hat)

library(reshape2)
library(ggplot2)
colnames(beta.hat) <- paste("beta", 1:16, sep = ".")
beta.df <- melt(as.data.frame(beta.hat))
pdf(file = "pen-16d-adj.pdf", width = 9,height = 5)
ggplot(beta.df ,aes(x = variable,y = value))  + geom_boxplot()
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
fdp <- sum(beta.hat[,-c(7,8)] != 0) / length(beta.hat[,-c(7,8)])
print(fdp)
fnp <- sum(beta.hat[,c(7,8)] == 0) / length(beta.hat[,c(7,8)])
print(fnp)
colMeans(beta.hat)
sum(colMeans(beta.hat))
mean(pars.hat)
sd(pars.hat)

library(reshape2)
library(ggplot2)
colnames(beta.hat) <- paste("beta", 1:16, sep = ".")
beta.df <- melt(as.data.frame(beta.hat))
pdf(file = "al-16d-adj.pdf", width = 9,height = 5)
ggplot(beta.df ,aes(x = variable,y = value))  + geom_boxplot()
dev.off()

##################################################
## Parse 500 nnls fitting results of two adjecent delays, sd = 0.02
## Simulation script:
## nnls-16d-2dadj.R
## Commit: 471e1186e954baed9300525f4de0a657bf994793
## Thu Apr  3 10:59:04 CDT 2014
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
fdp <- sum(beta.hat[,-c(7,8)] != 0) / length(beta.hat[,-c(7,8)])
print(fdp)
fnp <- sum(beta.hat[,c(7,8)] == 0) / length(beta.hat[,c(7,8)])
print(fnp)
colMeans(beta.hat)
sum(colMeans(beta.hat))
mean(pars.hat)
sd(pars.hat)

library(reshape2)
library(ggplot2)
colnames(beta.hat) <- paste("beta", 1:16, sep = ".")
beta.df <- melt(as.data.frame(beta.hat))
pdf(file = "nnls-16d-adj-sd02.pdf", width = 9,height = 5)
ggplot(beta.df ,aes(x = variable,y = value))  + geom_boxplot()
dev.off()
rm(list = ls())

##################################################
## Parse 500 penalized fitting results of two adjecent delays, sd = 0.02
## Simulation script:
## pen-16d-2dadj.R
## Commit: 45f5ae3a8567e5398b8b3048d85035d607b87460
## Thu Apr 10 09:32:03 CDT 2014
##################################################

pars.hat <- beta.hat <- c()
beta.true <- rep(0, 16)
beta.true[c(7,8)] <- 2
pars.true <- 0.5

for(i in 0:19){
    load(paste("pen-2dadj-sd02-", i, ".RData", sep=""))
    for(j in 1:25){
        pars.hat <- c(pars.hat, pen.res[[j]]$pars)
        beta.hat <- rbind(beta.hat, pen.res[[j]]$beta)
    }
}

cat("pen\n")
fdp <- sum(beta.hat[,-c(7,8)] != 0) / length(beta.hat[,-c(7,8)])
print(fdp)
fnp <- sum(beta.hat[,c(7,8)] == 0) / length(beta.hat[,c(7,8)])
print(fnp)
colMeans(beta.hat)
sum(colMeans(beta.hat))
mean(pars.hat)
sd(pars.hat)

library(reshape2)
library(ggplot2)
colnames(beta.hat) <- paste("beta", 1:16, sep = ".")
beta.df <- melt(as.data.frame(beta.hat))
pdf(file = "pen-16d-adj-sd02.pdf", width = 9,height = 5)
ggplot(beta.df ,aes(x = variable,y = value))  + geom_boxplot()
dev.off()
rm(list = ls())

##################################################
## Parse 500 adaptive-lasso fitting results of two adjecent delays, sd = 0.02
## Simulation script:
## al-16d-2dadj.R
## Commit: 45f5ae3a8567e5398b8b3048d85035d607b87460
## Thu Apr 10 09:32:03 CDT 2014
##################################################

pars.hat <- beta.hat <- c()
beta.true <- rep(0, 16)
beta.true[c(7,8)] <- 2
pars.true <- 0.5

for(i in 0:19){
    load(paste("al-2dadj-sd02-", i, ".RData", sep=""))
    for(j in 1:25){
        pars.hat <- c(pars.hat, al.res[[j]]$pars)
        beta.hat <- rbind(beta.hat, al.res[[j]]$beta)
    }
}

cat("Adaptive LASSO\n")
fdp <- sum(beta.hat[,-c(7,8)] != 0) / length(beta.hat[,-c(7,8)])
print(fdp)
fnp <- sum(beta.hat[,c(7,8)] == 0) / length(beta.hat[,c(7,8)])
print(fnp)
colMeans(beta.hat)
sum(colMeans(beta.hat))
mean(pars.hat)
sd(pars.hat)

library(reshape2)
library(ggplot2)
colnames(beta.hat) <- paste("beta", 1:16, sep = ".")
beta.df <- melt(as.data.frame(beta.hat))
pdf(file = "al-16d-adj-sd02.pdf", width = 9,height = 5)
ggplot(beta.df ,aes(x = variable,y = value))  + geom_boxplot()
dev.off()
rm(list = ls())


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
        g <- cbind(nnls.res[[j]]$Xdf, nnls.res[[j]]$Zdf)
        H <- t(g)%*%g
        count <- c
        Covar[[i*25 + j]] <- NeweyWest.Var( 0.5*(t(H)+H) ,g,10)
        allpars <- c(pars.hat, beta.hat)
        cover <- ((allpars + 1.96 * sqrt(diag(Covar[[i*25 + j]]))) >= c(pars.true, beta.true)) & ((allpars - 1.96 * sqrt(diag(Covar[[i*25 + j]]))) <= c(pars.true, beta.true))
        coverage <- rbind(coverage, cover)
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
pdf(file = "pen-16d-adj-sd02.pdf", width = 9,height = 5)
ggplot(beta.df ,aes(x = variable,y = value))  + geom_boxplot()
dev.off()
rm(list = ls())

load(paste("blowfly-nnls-500-1.RData", sep=""))
load("data-blowfly-1000.RData")

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


for(i in 0:19){
    load(paste("blowfly-nnls-500-", i, ".RData", sep=""))
    for(j in 1:25){
        pars.hat <- rbind(pars.hat, nnls.res[[j]]$pars)
        beta.hat <- rbind(beta.hat, nnls.res[[j]]$beta)
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

