source("./R/timevar.R")
## source("./sources/LS.sparse.R")
## source("./sources/DSIRfnSparse.R")
## source("./sources/poslasso.R")

library(penalized)
library(CollocInfer)
library(deSolve)
library(nnls)

load("sim.tv02.sd200.RData")
nnlist <- sim.res

times <- seq(0, 5, by = 1/52)
knots <- times
norder = 3
nbasis = length(knots) + norder - 2
range  = range(knots)
basis <- create.bspline.basis(range=range(knots), nbasis=nbasis, norder=norder, breaks=knots)

set.seed(42)
sim.res <- list()
for(i in 1:100){
    print(i)
    initPars <- nnlist[[i]]$res$pars
    initKappa <- nnlist[[i]]$res$kappa
    coefs <- nnlist[[i]]$res$coefs
    tv.fused <- LS.tv(tvDSIRfn, nnlist[[i]]$data, times, pars = initPars, kappa = initKappa, coefs = coefs, basisvals = basis, lambda = 1000, in.meth='nlminb', control.out = list(method = "fused", maxIter = 10, lambda.sparse = -1), nnls.res = nnlist[[i]]$res)
    sim.res[[i]] <- tv.fused$select
    save(sim.res, file ="fused-tv02-sd200.RData")
}
