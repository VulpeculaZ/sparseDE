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
