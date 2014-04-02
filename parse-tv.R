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

library(reshape2)
library(ggplot2)
colnames(kappa.nnls) <- paste("k", c(1:12), sep="")
kappa.df <- melt(as.data.frame(kappa.nnls))
pdf("tv-nnls-02.pdf", 7,4)
ggplot(kappa.df ,aes(x = variable,y = value))  + geom_boxplot()
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

## Including all
kappamat <- matrix(NA, 100, 12)
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
pdf("tv-fused-02.pdf", 7,4)
ggplot(kappa.df ,aes(x = variable,y = value))  + geom_boxplot() + ylim(0,0.02 )
dev.off()


##################################################
load("fused-tv02.RData")
