library(reshape2)
library(ggplot2)
colnames(kappa.fused1) <- paste("k", c(1:12), sep="")
kappa.df <- melt(as.data.frame(kappa.fused1))
pdf("tv-fused1-02.pdf", 9, 5)
ggplot(kappa.df ,aes(x = variable,y = value))  + geom_boxplot() + ylim(0,0.02 )
dev.off()



##################################################
## Parse 500 fused-lasso fitting results of 1 delays, 1 jump in tv, sd = 100
## Simulation script: delay-tv-lasso.R
## Fri Sep  5 11:28:00 CDT 2014
## Commit: 304a4e007e1f94556c7640826a117cc8175c3198
##################################################

pars.pen <- beta.pen <- kappa.pen <- c()
beta.true <- rep(0, 6)
beta.true[5] <- 1
pars.true <- 10
kappa.true <- c(rep(0.005,3), rep(0.0025,3))

for(i in 1:19){
    load(paste("tv-lasso-6d-6tv-sd100-", i, ".RData", sep=""))
    for(j in 1:25){
        pars.pen <- c(pars.pen, sim.res[[j]]$select$pars.pen)
        beta.pen <- rbind(beta.pen, sim.res[[j]]$select$beta.pen)
        kappa.pen <- rbind(kappa.pen, sim.res[[j]]$select$kappa.pen)
    }
}

cat("lasso and fused lasso")
fdp <- sum(beta.pen[,-c(6)] != 0) / length(beta.pen[,-c(5)])
print(fdp)
fnp <- sum(beta.pen[,c(6)] == 0) / length(beta.pen[,c(5)])
print(fnp)
colMeans(beta.pen)

fdp <- sum((kappa.pen[,-6] - kappa.pen[,-1])[,-3] != 0 ) / length((kappa.pen[,-6] - kappa.pen[,-1])[,-3])
fdp
fnp <- sum(kappa.pen[,3]-kappa.pen[,4] == 0) / dim(kappa.pen)[1]
fnp
colMeans(kappa.pen)

library(reshape2)
library(ggplot2)
colnames(beta.pen) <- paste("beta", 1:6, sep = ".")
beta.df <- melt(as.data.frame(beta.pen))
pdf(file = "lasso-6d-6tv-beta.pdf", width = 9,height = 5)
ggplot(beta.df ,aes(x = variable,y = value))  + geom_boxplot()
dev.off()



colnames(kappa.pen) <- paste("kappa", 1:6, sep = ".")
kappa.df <- melt(as.data.frame(kappa.pen))
kappa.df$value[kappa.df$value > 0.02 | kappa.df$value < 0.001]  <- NA
pdf(file = "lasso-6d-6tv-kappa4.pdf", width = 9,height = 5)
ggplot(kappa.df ,aes(x = variable,y = value))  + geom_boxplot()
dev.off()

##################################################
## Parse 500 nnls fitting results of 1 delays, 1 jump in tv, sd = 100
## Simulation script: nnls-6d-6tv.R
## Script:
## nnls-6d-6tv.R
## Commit: 304a4e007e1f94556c7640826a117cc8175c3198
##################################################

pars.hat <- beta.hat <- kappa.hat <- c()
beta.true <- rep(0, 6)
beta.true[5] <- 1
pars.true <- 10
kappa.true <- c(rep(0.005,3), rep(0.0025,3))

for(i in 0:19){
    load(paste("nnls-6d-6tv-sd100-", i, ".RData", sep=""))
    for(j in 1:25){
        pars.hat <- c(pars.hat, nnls.res[[j]]$pars)
        beta.hat <- rbind(beta.hat, nnls.res[[j]]$conv$pars.kappa.beta[dim(nnls.res[[j]]$conv$pars.kappa.beta)[1], 8:13]
)
        kappa.hat <- rbind(kappa.hat, nnls.res[[j]]$kappa)
    }
}

cat("nnls")
fdp <- sum(beta.hat[,-c(6)] != 0) / length(beta.hat[,-c(6)])
print(fdp)
fnp <- sum(beta.hat[,c(6)] == 0) / length(beta.hat[,c(6)])
print(fnp)
colMeans(beta.hat)


fdp <- sum((kappa.hat[,-6] - kappa.hat[,-1])[,-3] != 0 ) / length((kappa.hat[,-6] - kappa.hat[,-1])[,-3])
fdp
fnp <- sum(kappa.hat[,3]-kappa.hat[,4] == 0) / dim(kappa.hat)[1]
fnp

colMeans(kappa.hat)

library(reshape2)
library(ggplot2)
colnames(beta.hat) <- paste("beta", 1:6, sep = ".")
beta.df <- melt(as.data.frame(beta.hat))
pdf(file = "nnls-6d-6tv-beta.pdf", width = 9,height = 5)
ggplot(beta.df ,aes(x = variable,y = value))  + geom_boxplot()
dev.off()

colnames(kappa.hat) <- paste("kappa", 1:6, sep = ".")
kappa.df <- melt(as.data.frame(kappa.hat))
pdf(file = "nnls-6d-6tv-kappa.pdf", width = 9,height = 5)
ggplot(kappa.df, aes(x = variable,y = value))  + geom_boxplot()
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

