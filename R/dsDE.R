##################################################
## Drift Spline for Differential Equations
##################################################

dsde.LS <- function(fn, data, times, pars, coefs = NULL, basisvals = NULL,
    lambda, fd.obj = NULL, more = NULL, weights = NULL, quadrature = NULL,
    in.meth = "nlminb", control.in = list(),
    eps = 1e-06, active = NULL, posproc = FALSE,
    poslik = FALSE, discrete = FALSE, names = NULL, sparse = FALSE,
    likfn = make.id(), likmore = NULL)
{
    if (is.null(active)) {
        active = 1:length(pars)
    }
    profile.obj <- LS.setup(pars = pars, coefs = coefs, fn = fn,
                            basisvals, lambda = lambda, fd.obj, more, data, weights,
                            times, quadrature, eps = 1e-06, posproc, poslik, discrete,
                            names, sparse, likfn = make.id(), likmore = NULL)
    dims = dim(data)
    lik = profile.obj$lik
    proc = profile.obj$proc
    proc$d2fdp2 <- d2fdp2
    coefs = profile.obj$coefs
    data = profile.obj$data
    times = profile.obj$times
    coefs.pars <- c(coefs, pars)
    ## Optimization:
        check.lik.proc.data.coefs(lik, proc, data, times, coefs)
    if (in.meth == "optim") {
        if (is.null(control.in$trace)) {
            control.in$trace = 0
        }
        if (is.null(control.in$maxit)) {
            control.in$maxit = 1000
        }
        if (is.null(control.in$reltol)){
            control.in$reltol = 1e-12
        }
        if (is.null(control.in$meth)) {
            control.in$meth = "BFGS"
        }
        imeth = control.in$meth
        control.in$meth = NULL
        res = optim(coefs.pars, SplineCoefsErrAll, gr = SplineCoefsDAll,
        hessian = TRUE, control = control.in,
        times = times, data = data, lik = lik, proc = proc,
        method = imeth)
        ncoefs = matrix(res$par[1:length(coefs)], ncol(lik$bvals), length(res$par)/ncol(lik$bvals))
        pars <- res$par[(1+length(coefs)) : length(res$par)]
    }
    else if (in.meth == "nlminb") {
        if (is.null(control.in$trace)) {
            control.in$trace = 0
        }
        if (is.null(control.in$eval.max)) {
            control.in$eval.max = 2000
        }
        if (is.null(control.in$iter.max)) {
            control.in$iter.max = 1000
        }
        if (is.null(control.in$rel.tol)) {
            control.in$rel.tol = 1e-12
        }
        res <- nlminb(coefs.pars, SplineCoefsErrAll, gradient = SplineCoefsDAll,
                      hessian = LS.Hessian, control = control.in, times = times,
                      data = data, lik = lik, proc = proc)
        ncoefs <- matrix(res$par[1:length(coefs)], ncol(lik$bvals), length(res$par)/ncol(lik$bvals))
        pars <- res$par[(1+length(coefs)) : length(res$par)]
        names(pars) <- proc$more$parnames
    }
    else {
        stop("Unknown optimizer specified")
    }
    if (!is.null(proc$more$names)) {
        colnames(ncoefs) = proc$more$names
    }
    return(list(coefs = ncoefs, pars = pars, res = res, data = data))
}

SplineCoefsErrAll <- function(coefs.pars, times, data, lik, proc, sgn = 1){
    coefs <- coefs.pars[1:(length(coefs.pars) - length(proc$more$parnames))]
    pars <- coefs.pars[(length(coefs.pars) - length(proc$more$parnames) + 1):length(coefs.pars)]
    names(pars) <- proc$more$parnames
    res <- SplineCoefsErr(coefs, times, data, lik, proc, pars, sgn)
    return(res)
}

SplineCoefsDAll <- function(coefs.pars, times, data, lik, proc,  sgn = 1){
    coefs <- coefs.pars[1:(length(coefs.pars) - length(proc$more$parnames))]
    pars <- coefs.pars[(length(coefs.pars) - length(proc$more$parnames) + 1):length(coefs.pars)]
    names(pars) <- proc$more$parnames
    DJDC <- SplineCoefsDC(coefs, times, data, lik, proc, pars, sgn)
    DJDP <- SplineCoefsDP(coefs, times,  data, lik, proc, pars, sgn)
    return(c(DJDC, DJDP))
}

LS.Hessian <- function(coefs.pars, times, data, lik, proc, pars, sgn = 1){
    coefs <- coefs.pars[1:(length(coefs.pars) - length(proc$more$parnames))]
    pars <- coefs.pars[(length(coefs.pars) - length(proc$more$parnames) + 1):length(coefs.pars)]
    names(pars) <- proc$more$parnames
    DJDC2 <- SplineCoefsDC2(coefs, times, data, lik, proc, pars, sgn)
    DJDCDP <- SplineCoefsDCDP(coefs, times, data, lik, proc, pars, sgn)
    DJDP2 <- SplineCoefsDP2(coefs, times, data, lik, proc, pars, sgn)

    H <- cbind(DJDC2, DJDCDP)
    H <- rbind(H, cbind(t(DJDCDP), DJDP2))
    return(H)
}



SplineCoefsDP2 <- function(coefs, times, data, lik, proc, pars,  sgn = 1){
    coefs2 = matrix(coefs, ncol(lik$bvals), length(coefs)/ncol(lik$bvals))
    H <- proc$d2fdp2(coefs2, proc$bvals, pars, proc$more)
    return(as.matrix(sgn * H))
}

d2fdp2 <- function (coefs, bvals, pars, more){
    devals = as.matrix(bvals$bvals %*% coefs)
    ddevals = as.matrix(bvals$dbvals %*% coefs)
    colnames(devals) = more$names
    colnames(ddevals) = more$names
    fdevals <- more$fn(more$qpts, devals, pars, more$more)
    difs <- ddevals - fdevals
    dfdp <- more$dfdp(more$qpts, devals, pars, more$more)
    d2fdp2 <- more$d2fdp2(more$qpts, devals, pars, more$more)
    weights <- checkweights(more$weights, more$whichobs,dfdp[,,1, drop = FALSE])
    H <- matrix(, nrow = length(pars), ncol = length(pars))
    for(i in 1:length(pars)){
        for(j in 1:length(pars)){
            H[i,j] <- sum((dfdp[,,i] * dfdp[,,j] - d2fdp2[,,i,j] * difs) * weights)
        }
    }
    return(2 * H)
}
