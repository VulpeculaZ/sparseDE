##################################################
## Parse 500 nnls fitting results of two adjecent delays
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

