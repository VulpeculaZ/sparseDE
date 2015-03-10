ProfileSSE.covariance.delay <- function(pars, beta, active = NULL, eps = 1e-06, ...)
{
    if (is.null(active)) {
        active = 1:length(pars)
    }
    apars <- pars[active]
    H <- matrix(0, length(apars) + length(beta), length(apars) + length(beta))
    g <- ProfileDP.sparse(pars = pars, beta = beta, active = active, ...)
    gg <- c(colSums(g$Xdf), colSums(g$Zdf))
    for(i in 1:(length(apars) + length(beta))){
        if(i <= length(apars)){
            tpars <- pars
            tpars[active][i] <-tpars[active][i] + eps
            tg <- ProfileDP.sparse(tpars, beta, active = active, ...)
            tg <- c(colSums(tg$Xdf), colSums(tg$Zdf))
        } else {
            tbeta <- beta
            tbeta[i - length(apars)] <- beta[i - length(apars)] + eps
            tbeta <- tbeta / sum(tbeta)
            tg <- ProfileDP.sparse(pars, tbeta, active = active, ...)
            tg <- c(colSums(tg$Xdf), colSums(tg$Zdf))
        }
        H[,i] <- (tg - gg)/eps
    }
    Covar <- NeweyWest.Var( 0.5*(t(H)+H) ,cbind(g$Xdf, g$Zdf) ,5)
    return(Covar)
}


ProfileDP.sparse <- function(pars, beta, fn, data, times, coefs = NULL, basisvals = NULL,
    lambda, fd.obj = NULL, more = NULL, weights = NULL, quadrature = NULL,
    in.meth = "nlminb", out.meth = "nls", control.in = list(),
    control.out = list(), eps = 1e-06, active = NULL, posproc = FALSE,
    poslik = FALSE, discrete = FALSE, names = NULL, sparse = FALSE,
    likfn = make.id(), likmore = NULL, delay = NULL, tauMax = NULL,
    basisvals0 = NULL, coefs0 = NULL, nbeta, ndelay, tau)
{
    allpars <- pars
    betanames <- c()
    for(i in 1:length(nbeta)){
        for(j in 1:nbeta[i]){
            betanames <- c(betanames,paste("beta",i,".",j, sep = ""))
        }
    }
    ## Create y.d
    fdnames <- list(NULL, NULL, NULL)
    fdnames[[2]] <- attr(coefs, "dimnames")[[2]]
    fdobj0 <- list(coefs = coefs0, basis = basisvals0, fdnames =fdnames)
    fdobj.d <- list(coefs = coefs, basis = basisvals, fdnames =fdnames)
    attr(fdobj0, "class") <- "fd"
    attr(fdobj.d, "class") <- "fd"
    profile.obj = LS.setup(pars = pars, coefs = coefs, fn = fn,
    basisvals, lambda = lambda, fd.obj, more, data, weights,
        times, quadrature, eps = 1e-06, posproc, poslik, discrete,
        names, sparse, likfn = make.id(), likmore = NULL)
    dims = dim(data)
    lik = profile.obj$lik
    proc = profile.obj$proc
    coefs = profile.obj$coefs
    data = profile.obj$data
    times = profile.obj$times
    proc$more$betanames <- betanames
    ##################################################
    ## Added delay data and functions
    ##################################################
    delayProcObj <- delay.fit.sparse(fd0 = fdobj0, fd.d = fdobj.d, times = proc$more$qpts, tau = tau, beta= beta, ndelay = ndelay, in.meth = in.meth)
    delayLikObj <- delay.fit.sparse(fd0 = fdobj0, fd.d = fdobj.d, times = times,tau = tau, beta= beta, ndelay = ndelay, in.meth = in.meth)
    lik$more$more$y.d <- delayLikObj$y.d
    proc$more$more$y.d <- delayProcObj$y.d
    lik$more$more$bvals.d <- delayLikObj$bvals.d
    proc$more$more$bvals.d <- delayProcObj$bvals.d
    proc$more$more$bvals.d.list <- delayProcObj$bvals.d.list
    proc$more$more$y.d.list <- delayProcObj$y.d.list
    proc$more$more$ndelay <- lik$more$more$ndelay <- ndelay
    proc$more$more$nbeta <- lik$more$more$nbeta <- sapply(tau, length)
    proc$more$more$tau <- lik$more$more$tau <- tau
    proc$dfdc <- delay$dfdc
    proc$d2fdc2 <- delay$d2fdc2.DDE
    proc$d2fdcdp <- delay$d2fdcdp.sparse
    proc$more$delay <- delay
    proc$more$dfdtau <- dfdbeta.sparse
    proc$more$d2fdxdtau <- d2fxdbeta.sparse
    proc$more$d2fdx.ddtau <- d2fdx.ddbeta.sparse
    ## Ires <- IresTmp

    apars = pars[active]
    aparamnames = names(apars)
    if (is.null(control.out$maxIter)) {
        control.out$maxIter = 100
    }
    if (is.null(control.out$tol)){
        control.out$tol = 1e-08
    }
    control.in$iter.max = 1

    res <- ProfileDP.AllPar.sparse(pars = pars, beta = beta, times = times, data = data, coefs = coefs, lik = lik, proc = proc, in.meth = in.meth, control.in = control.in, basisvals = basisvals, fdobj0 = fdobj0)
    Xdf <- -res$Xdf
    Zdf <- -res$Zdf
    return(list(Xdf = Xdf, Zdf = Zdf))
}

ProfileDP.AllPar.sparse <- function(pars, beta, times, data, coefs, lik, proc,
                             in.meth='nlminb', control.in=NULL,
                             dcdp=NULL, oldpars=NULL, use.nls=TRUE, sgn=1,
                             basisvals, fdobj0)
{
    fdnames <- list(NULL, NULL, NULL)
    fdnames[[2]] <- attr(coefs, "dimnames")[[2]]
    fdobj.d <- list(coefs = coefs, basis = basisvals, fdnames =fdnames)
    attr(fdobj.d, "class") <- "fd"
    delayProcObj <- delay.fit.sparse(fd0 = fdobj0, fd.d = fdobj.d, times = proc$more$qpts, tau = proc$more$more$tau, beta= beta, ndelay = proc$more$more$ndelay, in.meth = in.meth)
    delayLikObj <- delay.fit.sparse(fd0 = fdobj0, fd.d = fdobj.d, times = times,tau = lik$more$more$tau, beta= beta, ndelay = lik$more$more$ndelay, in.meth = in.meth)
    lik$more$more$y.d <- delayLikObj$y.d
    proc$more$more$y.d <- delayProcObj$y.d
    lik$more$more$bvals.d <- delayLikObj$bvals.d
    proc$more$more$bvals.d <- delayProcObj$bvals.d
    proc$more$more$bvals.d.list <- delayProcObj$bvals.d.list
    proc$more$more$y.d.list <- delayProcObj$y.d.list
    ## Calculate fitted value after inner optimization:
    devals = as.matrix(lik$bvals%*%coefs)
    colnames(devals) = proc$more$names
    ## Squared errors: No need to change for DDE
    weights = checkweights(lik$more$weights,lik$more$whichobs,data)
    dlikdp = lik$more$dfdp(times,devals,pars,lik$more$more)
    dlikdp = matrix(dlikdp,dim(dlikdp)[1]*dim(dlikdp)[2],dim(dlikdp)[3])
    ## dlikdp will be zero if the likelihood doesn't directly have ode parameters,
    ## which is true for least square case.
    dlikdx = lik$more$dfdx(times,devals,pars,lik$more$more)
    ## dlikdx[i,,] is an identity matrix for every i.
    dlikdc = c()
    for(i in 1:dim(dlikdx)[2]){
        tH = c()
        for(j in 1:dim(dlikdx)[3]){
            ## Only dlikdx[,i,i] are non-zero (all 1's)
            tH = cbind(tH,as.matrix(diag(dlikdx[,i,j])%*%lik$bvals))
        }
        dlikdc = rbind(dlikdc,tH)
    }
    ## ??dlikdc: why 0.5 and 0.5 ??
    d2Hdc2  = SplineCoefsDC2sparse(coefs,times,data,lik,proc,pars)
    ## need not be changed?
    d2Hdcdp = SplineCoefsDCDP.sparse(coefs, times, data, lik, proc, pars)
    ## Got warning message:
    ## In dim(weights[whichrows, ]) == dim(diffs):
    ## longer object length is not a multiple of shorter object length?

    ## Use Implicit function theorem:
    if(is.matrix(d2Hdc2)){
        ## When will it not be a matrix? How to use solve in that case?
        dcdp = ginv(d2Hdc2) %*% d2Hdcdp
    } else {
        dcdp = as.matrix(solve(d2Hdc2,d2Hdcdp))
    }
    ## Chain rule:
    df = dlikdc%*%dcdp ## + dlikdp
    Xdf <- df[, 1:length(pars), drop = FALSE]
    Zdf <- df[, (length(pars) + 1): dim(df)[2]]
    return(list(Xdf = Xdf, Zdf = Zdf))
}
