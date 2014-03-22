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
