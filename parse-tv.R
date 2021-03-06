load("sim.tv01.RData")
tv.nnls <- sim.res
rm(list = "sim.res")
tv.fused <- sim.res
rm(list = "sim.res")

pars.true <- 0.5
kappa.true <- rep(0.01, 12)
kappa.true[10:12] <- 0.006
kappa.nnls <- gamma.nnls <- kappa.fused <- gamma.fused <- c()

for(i in 1:length(tv.nnls)){
    kappa.nnls <- rbind(kappa.nnls, tv.nnls[[i]]$res$kappa)
    gamma.nnls <- c(gamma.nnls, tv.nnls[[i]]$res$pars)
    kappa.fused <- rbind(kappa.fused, tv.fused[[i]]$kappa)
    gamma.fused <- c(gamma.fused, tv.fused[[i]]$pars)
}

pdf("barnnls.pdf", 7,6)
tmp <- colMeans(kappa.nnls)
names(tmp) <- names(colMeans(kappa.fused))
barplot(tmp)
dev.off()
pdf("barfused.pdf", 7,6)
barplot(colMeans(kappa.fused))
dev.off()
mean(sapply(apply(kappa.fused, 1, unique), length)) / 12


pars.pen <- kappa.pen <- coefs.pen <- list()

for(i in 1:length(lambda)){
    lambda.sparse <- lambda[i]
    res.sparse <- penalized(response = tv.nnls[[i]]$res$y, penalized = tv.nnls[[i]]$res$Zdf, unpenalized = tv.nnls[[i]]$res$Xdf, lambda1=0, lambda2 = lambda[i], positive = TRUE, fusedl=TRUE)
    pars.pen[[i]] <- res.sparse@unpenalized
    kappa.pen[[i]] <- res.sparse@penalized
}

kappa.pen

colSums(kappa.nnls)
apply(kappa.nnls, 2,sd)

tmp <- sweep(tv.nnls[[1]]$res$Zdf, 1, mean(tv.nnls[[1]]$res$Zdf))


##################################################
## Parse 100 nnls fitting results of estimating time varying coefficients
## Simulation script:
## sim.tv.nnls.R
## Commit: df114f8982bc1c11c228a4de5ca507198df08e66
## Tue Apr  1 13:36:23 CDT 2014
##################################################

load("sim.tv02.RData")
tv.nnls <- sim.res
pars.true <- 10
kappa.true <- rep(0.005, 12)
kappa.true[10:12] <- 0.002

kappa.nnls <- gamma.nnls  <- c()
for(i in 1:length(tv.nnls)){
    kappa.nnls <- rbind(kappa.nnls, tv.nnls[[i]]$res$kappa)
    gamma.nnls <- c(gamma.nnls, tv.nnls[[i]]$res$pars)
}

mean(gamma.nnls)
sd(gamma.nnls)
library(reshape2)
library(ggplot2)
colnames(kappa.nnls) <- paste("k", c(1:12), sep="")
kappa.df <- melt(as.data.frame(kappa.nnls))
pdf("tv-nnls-02.pdf", 9,5)
ggplot(kappa.df ,aes(x = variable,y = value))  + geom_boxplot() + ylim(0,0.02)
dev.off()



##################################################
## Parse 100 fused-lasso fitting results of estimating time varying coefficients
## Fused-lasso was run at every outer-iteration
## Simulation script:
## sim.tv.fused02.R
## Commit: df114f8982bc1c11c228a4de5ca507198df08e66
## Tue Apr  1 14:23:24 CDT 2014
##################################################
library(CollocInfer)

pars.hat <- kappa.hat <- c()
kappa.true <- rep(0.005, 12)
kappa.true[c(10,11,12)] <- 0.002
pars.true <- 10

times <- seq(0, 5, by = 1/52)
knots <- times
norder = 3
nbasis = length(knots) + norder - 2
range  = range(knots)
basis <- create.bspline.basis(range=range(knots), nbasis=nbasis, norder=norder, breaks=knots)
bvals <- eval.basis(times, basis, 0)
bic <- matrix(NA, 50, 100)


for(i in 1:50){
    if(!file.exists(paste("lassotv", i, ".RData", sep=""))){
        cat("Iteration:",i, "does not exist", "\n")
    }else{
        load(paste("lassotv", i, ".RData", sep=""))
        if(length(sim.res) != 100)
        cat("Iteration prematurally terminated:",i, "\n")
        for(j in 1:length(sim.res)){
            devals <- as.matrix(bvals%*%sim.res[[j]]$ncoefs)
            f <- sim.res[[j]]$data - devals
            sd.pen <- sd(f)
            ll.pen <- - sum(f^2) / (sd.pen^2) / 2 - length(f) * log(sd.pen)
            bic[i,j] <- -2 * ll.pen + (length(unique(sim.res[[j]]$res$kappa)) + length(sim.res[[j]]$res$pars)) * log(length(f))
        }
    }
}
minbic <- apply(bic, 2, min, na.rm = "TRUE")
posminbic <- apply(sweep(bic, 2 , minbic, "=="),2,which)

## Including all lambda's
kappaMat <- matrix(NA, 100, 12)
pars <- c()
fdp.kappa <- fnp.kappa <- 0
for(i in 1:100){
    filenamei <- paste("lassotv", posminbic[i],".RData", sep = "")
    load(filenamei)
    kappaMat[i,] <- sim.res[[i]]$res$kappa
    pars <- c(pars, sim.res[[i]]$res$pars)
}

fdp <- sum((kappaMat[,-12] - kappaMat[,-1])[,-9] != 0 ) / length((kappaMat[,-12] - kappaMat[,-1])[,-9])
fnp <- sum(kappaMat[,9]-kappaMat[,10] == 0) / dim(kappaMat)[1]
mean(pars)

## Try plotting at log scale
library(reshape2)
library(ggplot2)
colnames(kappaMat) <- paste("k", c(1:12), sep="")
kappa.df <- melt(as.data.frame(kappaMat))
pdf("tv-fused-02.pdf", 9,5)
ggplot(kappa.df ,aes(x = variable,y = value))  + geom_boxplot() + ylim(0,0.02 )
dev.off()


##################################################
## Parse 100 fused-lasso fitting results of estimating time varying coefficients
## Fused-lasso was run at outer-iteration nnls converged
## Simulation script:
## fused-tv-02.R
## Wed Apr  2 14:18:36 CDT 2014
##################################################
load("fused-tv02.RData")
tv.fused1 <- sim.res
pars.true <- 10
kappa.true <- rep(0.005, 12)
kappa.true[10:12] <- 0.002

kappa.fused1 <- gamma.fused1  <- c()
for(i in 1:length(tv.fused1)){
    kappa.fused1 <- rbind(kappa.fused1, tv.fused1[[i]]$kappa)
    gamma.fused1 <- c(gamma.fused1, tv.fused1[[i]]$pars)
}

fdp <- sum((kappa.fused1[,-12] - kappa.fused1[,-1])[,-9] != 0 ) / length((kappa.fused1[,-12] - kappa.fused1[,-1])[,-9])
fnp <- sum(kappa.fused1[,9]-kappa.fused1[,10] == 0) / dim(kappa.fused1)[1]
mean(gamma.fused1)
sd(gamma.fused1)

library(reshape2)
library(ggplot2)
colnames(kappa.fused1) <- paste("k", c(1:12), sep="")
kappa.df <- melt(as.data.frame(kappa.fused1))
pdf("tv-fused1-02.pdf", 9, 5)
ggplot(kappa.df ,aes(x = variable,y = value))  + geom_boxplot() + ylim(0,0.02 )
dev.off()


##################################################
## Parse 100 nnls fitting results of estimating time varying coefficients
## Noise sd = 100
## Simulation script:
## sim.nnls.tv.sd100.R
## Commit: 471e1186e954baed9300525f4de0a657bf994793
## Thu Apr 10 10:22:42 CDT 2014
##################################################

load("sim.tv02.sd100.RData")
tv.nnls <- sim.res
pars.true <- 10
kappa.true <- rep(0.005, 12)
kappa.true[10:12] <- 0.002

kappa.nnls <- gamma.nnls  <- c()
for(i in 1:length(tv.nnls)){
    kappa.nnls <- rbind(kappa.nnls, tv.nnls[[i]]$res$kappa)
    gamma.nnls <- c(gamma.nnls, tv.nnls[[i]]$res$pars)
}
fdp <- sum((kappa.nnls[,-12] - kappa.nnls[,-1])[,-9] != 0 ) / length((kappa.nnls[,-12] - kappa.nnls[,-1])[,-9])
fdp
fnp <- sum(kappa.nnls[,9]-kappa.nnls[,10] == 0) / dim(kappa.nnls)[1]
fnp
mean(gamma.nnls)
sd(gamma.nnls)

library(reshape2)
library(ggplot2)
colnames(kappa.nnls) <- paste("k", c(1:12), sep="")
kappa.df <- melt(as.data.frame(kappa.nnls))
pdf("tv-nnls-sd100.pdf", 9,5)
ggplot(kappa.df ,aes(x = variable,y = value))  + geom_boxplot()
dev.off()


##################################################
## Parse 100 nnls fitting results of estimating time varying coefficients
## Noise sd = 200
## Simulation script:
## sim.nnls.tv.sd200.R
## Commit: 471e1186e954baed9300525f4de0a657bf994793
## Thu Apr 10 10:22:31 CDT 2014
##################################################

load("sim.tv02.sd200.RData")
tv.nnls <- sim.res
pars.true <- 10
kappa.true <- rep(0.005, 12)
kappa.true[10:12] <- 0.002

kappa.nnls <- gamma.nnls  <- c()
for(i in 1:length(tv.nnls)){
    kappa.nnls <- rbind(kappa.nnls, tv.nnls[[i]]$res$kappa)
    gamma.nnls <- c(gamma.nnls, tv.nnls[[i]]$res$pars)
}
fdp <- sum((kappa.nnls[,-12] - kappa.nnls[,-1])[,-9] != 0 ) / length((kappa.nnls[,-12] - kappa.nnls[,-1])[,-9])
fdp
fnp <- sum(kappa.nnls[,9]-kappa.nnls[,10] == 0) / dim(kappa.nnls)[1]
fnp
mean(gamma.nnls)
sd(gamma.nnls)

library(reshape2)
library(ggplot2)
colnames(kappa.nnls) <- paste("k", c(1:12), sep="")
kappa.df <- melt(as.data.frame(kappa.nnls))
pdf("tv-nnls-sd200.pdf", 9,5)
ggplot(kappa.df ,aes(x = variable,y = value))  + geom_boxplot() #+ ylim(0,0.01)
dev.off()


##################################################
## Parse 100 fused-lasso fitting results of estimating time varying coefficients
## Noise sd = 100
## Fused-lasso was run at outer-iteration nnls converged
## Simulation script:
## fused-tv-02.R
## Commit: 45f5ae3a8567e5398b8b3048d85035d607b87460
## Sun Apr 13 20:33:43 CDT 2014
##################################################
load("fused-tv02-sd100.RData")
tv.fused1 <- sim.res
pars.true <- 10
kappa.true <- rep(0.005, 12)
kappa.true[10:12] <- 0.002

kappa.fused1 <- gamma.fused1  <- c()
for(i in 1:length(tv.fused1)){
    kappa.fused1 <- rbind(kappa.fused1, tv.fused1[[i]]$kappa)
    gamma.fused1 <- c(gamma.fused1, tv.fused1[[i]]$pars)
}

fdp <- sum((kappa.fused1[,-12] - kappa.fused1[,-1])[,-9] != 0 ) / length((kappa.fused1[,-12] - kappa.fused1[,-1])[,-9])
fnp <- sum(kappa.fused1[,9]-kappa.fused1[,10] == 0) / dim(kappa.fused1)[1]
mean(gamma.fused1)
sd(gamma.fused1)

library(reshape2)
library(ggplot2)
colnames(kappa.fused1) <- paste("k", c(1:12), sep="")
kappa.df <- melt(as.data.frame(kappa.fused1))
pdf("tv-fused1-sd100.pdf", 9, 5)
ggplot(kappa.df ,aes(x = variable,y = value))  + geom_boxplot() + ylim(0,0.02 )
dev.off()


##################################################
## Parse 200 fused-lasso fitting results of estimating time varying coefficients
## Noise sd = 200
## Fused-lasso was run at outer-iteration nnls converged
## Simulation script:
## fused-tv-02.R
## Commit: 45f5ae3a8567e5398b8b3048d85035d607b87460
## Sun Apr 13 20:33:43 CDT 2014
##################################################
load("fused-tv02-sd200.RData")
tv.fused1 <- sim.res
pars.true <- 10
kappa.true <- rep(0.005, 12)
kappa.true[10:12] <- 0.002

kappa.fused1 <- gamma.fused1  <- c()
for(i in 1:length(tv.fused1)){
    kappa.fused1 <- rbind(kappa.fused1, tv.fused1[[i]]$kappa)
    gamma.fused1 <- c(gamma.fused1, tv.fused1[[i]]$pars)
}

fdp <- sum((kappa.fused1[,-12] - kappa.fused1[,-1])[,-9] != 0 ) / length((kappa.fused1[,-12] - kappa.fused1[,-1])[,-9])
fnp <- sum(kappa.fused1[,9]-kappa.fused1[,10] == 0) / dim(kappa.fused1)[1]
mean(gamma.fused1)

library(reshape2)
library(ggplot2)
colnames(kappa.fused1) <- paste("k", c(1:12), sep="")
kappa.df <- melt(as.data.frame(kappa.fused1))
pdf("tv-fused1-sd200.pdf", 9, 5)
ggplot(kappa.df ,aes(x = variable,y = value))  + geom_boxplot() + ylim(0,0.02 )
dev.off()
