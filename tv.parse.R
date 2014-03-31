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

## parse the result for nnls fit for time varying coefficients model:
## simulation done in: sim.tv.nnls.R
## Git SHA-1 fe560f681f791bb880fafb07b471dd714997f937
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

pdf("tmp.pdf", 7,6)
tmp <- colMeans(kappa.nnls)
## names(tmp) <- names(colMeans(kappa.fused))
barplot(tmp)
dev.off()

mean(gamma.nnls)


##################################################
## Parse 500 nnls fitting results of two adjecent delays
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
