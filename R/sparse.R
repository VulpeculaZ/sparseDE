nls.sparse <- function(pars, beta, active, basisvals, fdobj0, times, data, coefs, lik, proc, start, X.index, control.out, control.in, in.meth){
    if(control.out$method == "twoStage"){
        delta <- rep(1, length(data))
    }
    pars.names <- names(pars)
    f.conv <- pars.beta <- c()
    maxStep <- 8
    lambda.sparse <- control.out$lambda.sparse
    for(i in 1:control.out$maxIter){
        for(j in 1:maxStep){
            linObj <- ProfileSSE.AllPar.sparse(pars = pars, beta = beta, times = times, data = data, coefs = coefs, lik = lik, proc = proc, in.meth = in.meth,control.in = control.in, basisvals = basisvals, fdobj0 = fdobj0)
            f.new <- linObj$f
            f.new <- sum(f.new^2)
            if(control.out$method == "penalized"){
                f.new <- f.new + lambda.sparse * sum(abs(beta))
            }
            if(control.out$method == "enet"){
                f.new <- f.new + lambda.sparse * sum(abs(beta)) + 0.0005* sum(beta^2)
            }
            ## if(control.out$echo == TRUE){
            ## print(x = c(paste("Iter:", i, f.new)))
            ##    cat(pars, beta, "\n")
            ## }
            if(i == 1){
                break
            }else{
                if(f.conv[i - 1] - f.new > 0 & f.conv[i - 1] - f.new < control.out$tol){
                    return(list(pars=pars, beta = beta, coefs = coefs, f = f.new, y = y, Xdf = Xdf, Zdf = Zdf, conv = list(f = f.conv, pars.beta=pars.beta, conv.message = "Converged.")))
                }
                if(f.conv[i - 1] - f.new > 0){
                    break
                }
                if(f.conv[i - 1] - f.new < 0 & j == maxStep){
                    return(list(pars=pars.old, beta = beta.old, coefs = coefs, f = f.new,  y = y, Xdf = Xdf, Zdf = Zdf, conv = list(f = f.conv, pars.beta=pars.beta, conv.message = "Non-dereasing objective.")))
                }
                pars <- 0.5*(pars - pars.old) + pars.old
                beta <- 0.5*(beta - beta.old) + beta.old
            }
        }
        pars.old <- pars
        beta.old <- beta
        f.conv <- c(f.conv, f.new)
        pars.beta <- rbind(pars.beta, c(pars, beta))
        Xdf <- - linObj$df[, 1:length(pars), drop = FALSE]
        Zdf <- - linObj$df[, (length(pars) + 1): dim(linObj$df)[2]]
        y <- - linObj$df %*% c(pars, beta) + linObj$f
        coefs <- linObj$coefs
        if(control.out$method == "twoStage"){
            res <- twoStgEst(Z = Zdf, y = y, X = Xdf, delta = delta, n_lambda = 1, lambda0 = lambda.sparse)
            pars <- res$betamcp
            beta <- res$thetamcp
        }
        if(control.out$method == "penalized"){
            res <- penalized(response = y, penalized = Zdf, unpenalized = Xdf, lambda1 = lambda.sparse, positive = TRUE, trace = FALSE)
            pars <- res@unpenalized
            beta <- res@penalized
        }
        if(control.out$method == "enet"){
            res <- penalized(response = y, penalized = Zdf, unpenalized = Xdf, lambda1 = lambda.sparse, lambda2 = 0.00001, positive = TRUE, trace = FALSE)
            pars <- res@unpenalized
            beta <- res@penalized
        }
        if(control.out$method == "ols"){
            res <- lm.fit(x = cbind(Xdf, Zdf), y= y)
            pars <- res$coefficients[1:length(pars)]
            beta <- res$coefficients[(length(pars) + 1) : length(res$coefficients)]
        }
        if(control.out$method == "nnls"){
            res <- nnls(A = cbind(Xdf, Zdf), b= y)
            pars <- res$x[1:length(pars)]
            beta <- res$x[(length(pars) + 1) : length(res$x)]
        }
        names(pars) <- pars.names
    }
    return(list(pars=pars.old, beta = beta.old, coefs = coefs, f = f.new, y = y, Xdf = Xdf, Zdf = Zdf, conv = list(f = f.conv, pars.beta=pars.beta, conv.message = "Maximum iterations reached.")))
}



##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title Evaluate sparse delays.
##' @param fd0
##' @param fd.d
##' @param times
##' @param tau
##' @param beta
##' @param ndelay
##' @param basis
##' @return
##' @author Ziqian Zhou
delay.fit.sparse <- function(fd0, fd.d, times, tau, beta, ndelay, basis = NULL, lik = FALSE){
    basisvals0 <- fd0$basis
    basisvals.d <- fd.d$basis
    start.d <- fd.d$basis$rangeval[1]
    y.d.beta <- matrix(0, length(times), length(ndelay))
    j <- i <- 0
    y.d.list <- list()
    if(lik == TRUE){
        for(idelay in ndelay){
            j <- j + 1
            for(itau in tau[[j]]){
                i <- i + 1
                times.d <- times - itau
                y.d <- eval.fd(times.d[times.d >= start.d], fd.d)
                if(sum(times.d < start.d)){
                    y.d0 <- eval.fd(times.d[times.d < start.d], fd0)
                    y.d <- rbind(y.d0, y.d)
                }
                y.d.beta[,j] <- y.d.beta[,j] + beta[i] * y.d[,idelay]
            }
        }
        return(list(y.d = y.d.beta))
    }
    bvals.d <- bvals.d.list <- list()
    if(is.null(basis)){
        for(idelay in ndelay){
            j <- j + 1
            for(itau in tau[[j]]){
                i <- i + 1
                times.d <- times - itau
                y.d <- eval.fd(times.d[times.d >= start.d], fd.d)
                bvals <- eval.basis(times.d[times.d >= start.d], basisvals.d, 0)
                if(sum(times.d < start.d)){
                    y.d0 <- eval.fd(times.d[times.d < start.d], fd0)
                    y.d <- rbind(y.d0, y.d)
                    bvals <- rbind(matrix(0, nrow = length(times) - dim(bvals)[1], ncol = dim(bvals)[2]), bvals)
                }
                y.d.list[[i]] <- y.d[,idelay]
                bvals.d.list[[i]] <- bvals
                y.d.beta[,j] <- y.d.beta[,j] + beta[i] * y.d[,idelay]
                if(j == 1){
                    bvals.d[[j]] <- bvals.d.list[[i]] * beta[i]
                }
                else{
                    bvals.d[[j]] <- bvals.d[[j]] + bvals.d.list[[i]] * beta[i]
                }
            }
        }
    }
    else{
        bvals.d.list <- basis
        for(idelay in ndelay){
            j <- j + 1
            for(itau in tau[[j]]){
                i <- i + 1
                times.d <- times - itau
                y.d <- eval.fd(times.d[times.d >= start.d], fd.d)
                if(sum(times.d < start.d)){
                    y.d0 <- eval.fd(times.d[times.d < start.d], fd0)
                    y.d <- rbind(y.d0, y.d)
                }
                y.d.list[[i]] <- y.d[,idelay]
                y.d.beta[,j] <- y.d.beta[,j] + beta[i] * y.d[,idelay]
                if(j == 1){
                    bvals.d[[j]] <- bvals.d.list[[i]] * beta[i]
                }
                else{
                    bvals.d[[j]] <- bvals.d[[j]] + bvals.d.list[[i]] * beta[i]
                }
            }
        }
    }
    ## Names ??
    ## Returning only bvals.d[[1]] for now. Need to be fixed.
    return(list(y.d = y.d.beta, bvals.d = bvals.d[[1]], y.d.list = y.d.list, bvals.d.list = bvals.d.list))
}

inneropt.DDE <- function(data, times, pars, beta, coefs, lik, proc,
                         in.meth = "nlminb", control.in = list(),
                         basisvals, fdobj0)
{
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
        if (is.null(control.in$reportHessian)){
            control.in$reportHessian = TRUE
        }
        imeth = control.in$meth
        control.in$meth = NULL
        res = optim(coefs, SplineCoefsErr.DDE, gr = SplineCoefsDC.DDE,
        hessian = control.in$reportHessian, control = control.in,
        times = times, data = data, lik = lik, proc = proc, pars = pars,
        beta = beta, method = imeth, basisvals = basisvals, fdobj0 = fdobj0)
        ncoefs = matrix(res$par, ncol(lik$bvals), length(res$par)/ncol(lik$bvals))
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
        if (is.null(control.in$useHessian)) {
            Hessian = SplineCoefsDC2.DDE
        }
        else {
            Hessian = NULL
        }
        ## SplineCoefsErr do not need to be changed.
        ##
        res <- nlminb(coefs, SplineCoefsErr.DDE, gradient = SplineCoefsDC.DDE,
                      hessian = Hessian, control = control.in, times = times,
                      data = data, lik = lik, proc = proc, pars = pars,
                      basisvals = basisvals, fdobj0 = fdobj0, beta = beta)
        ncoefs = matrix(res$par, ncol(lik$bvals), length(res$par)/ncol(lik$bvals))
    }
    else {
        stop("Unknown optimizer specified")
    }
    if (!is.null(proc$more$names)) {
        colnames(ncoefs) = proc$more$names
    }
    return(list(coefs = ncoefs, res = res))
}

Profile.LS.sparse <- function(fn, data, times, pars, beta, coefs = NULL, basisvals = NULL,
    lambda, fd.obj = NULL, more = NULL, weights = NULL, quadrature = NULL,
    in.meth = "nlminb", out.meth = "nls", control.in = list(),
    control.out = list(), eps = 1e-06, active = NULL, posproc = FALSE,
    poslik = FALSE, discrete = FALSE, names = NULL, sparse = FALSE,
    likfn = make.id(), likmore = NULL, delay = NULL, tauMax = NULL,
    basisvals0 = NULL, coefs0 = NULL, nbeta, ndelay, tau)
{
    if (is.null(active)) {
        active = 1:length(pars)
    }
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
    delayProcObj <- delay.fit.sparse(fd0 = fdobj0, fd.d = fdobj.d, times = proc$more$qpts, tau = tau, beta= beta, ndelay = ndelay )
    delayLikObj <- delay.fit.sparse(fd0 = fdobj0, fd.d = fdobj.d, times = times,tau = tau, beta= beta, ndelay = ndelay)
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
    Ires <- inneropt.DDE(data, times, pars, beta = beta,  coefs, lik, proc, in.meth, control.in, basisvals = basisvals, fdobj0 = fdobj0)
    ## Ires <- IresTmp
    ncoefs <- Ires$coefs
    apars = pars[active]
    aparamnames = names(apars)
    if (is.null(control.out$maxIter)) {
        control.out$maxIter = 100
    }
    if (is.null(control.out$tol)){
        control.out$tol = 1e-08
    }
    res <- nls.sparse(pars = pars, beta = beta, active = active, basisvals = basisvals, fdobj0 = fdobj0, times = times, data = data, coefs = ncoefs, lik = lik, proc = proc, control.out = control.out, control.in = control.in, in.meth = in.meth)
    ncoefs <- res$coefs
    return(list( data = data,res = res, ncoefs = ncoefs))
}


ProfileSSE.AllPar.sparse <- function(pars, beta, times, data, coefs, lik, proc,
                             in.meth='nlminb', control.in=NULL,
                             dcdp=NULL, oldpars=NULL, use.nls=TRUE, sgn=1,
                             basisvals, fdobj0)
{
    ## Squared Error outer criterion
    ##    coefs = as.vector(coefs)  # First run the inner optimization
    ##
    ## f1 = SplineCoefsErr.DDE(coefs,times,data,lik,proc,pars, beta, basisvals = basisvals, fdobj0 = fdobj0)
    ## The coefficients for delayed times:
    ##################################################
    ## Not Sure:
    ## if(!is.null(dcdp)){
    ##     tcoefs = as.vector(coefs) + dcdp%*%(pars-oldpars);
    ##     f2 = SplineCoefsErr(tcoefs,times,data,lik,proc,pars)
    ##     if(f2 < f1){
    ##         coefs = tcoefs
    ##         f1 = f2
    ##     }
    ## }
    ###################################################
    ## Inner optimization need to be changed as well.
    Ires = inneropt.DDE(data,times,pars, beta, coefs,lik,proc, in.meth,control.in, basisvals = basisvals, fdobj0 = fdobj0)
    ncoefs = Ires$coefs
    fdnames <- list(NULL, NULL, NULL)
    fdnames[[2]] <- attr(coefs, "dimnames")[[2]]
    fdobj.d <- list(coefs = ncoefs, basis = basisvals, fdnames =fdnames)
    attr(fdobj.d, "class") <- "fd"
    ##################################################
    ## Added delay data and functions
    ##################################################
    delayProcObj <- delay.fit.sparse(fd0 = fdobj0, fd.d = fdobj.d, times = proc$more$qpts, tau = proc$more$more$tau, beta= beta, ndelay = proc$more$more$ndelay )
    delayLikObj <- delay.fit.sparse(fd0 = fdobj0, fd.d = fdobj.d, times = times,tau = lik$more$more$tau, beta= beta, ndelay = lik$more$more$ndelay)
    lik$more$more$y.d <- delayLikObj$y.d
    proc$more$more$y.d <- delayProcObj$y.d
    lik$more$more$bvals.d <- delayLikObj$bvals.d
    proc$more$more$bvals.d <- delayProcObj$bvals.d
    proc$more$more$bvals.d.list <- delayProcObj$bvals.d.list
    proc$more$more$y.d.list <- delayProcObj$y.d.list
    ## Calculate fitted value after inner optimization:
    devals = as.matrix(lik$bvals%*%ncoefs)
    colnames(devals) = proc$more$names
    ## Squared errors: No need to change for DDE
    weights = checkweights(lik$more$weights,lik$more$whichobs,data)
    f = as.vector(as.matrix(data - lik$more$fn(times, devals, pars, lik$more$more))*sqrt(weights))
    isnaf = is.na(f)
    f[isnaf] = 0
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
    d2Hdc2  = SplineCoefsDC2sparse(ncoefs,times,data,lik,proc,pars)
    ## need not be changed?
    d2Hdcdp = SplineCoefsDCDP.sparse(ncoefs, times, data, lik, proc, pars)
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
    df[isnaf,] = 0
    colnames(df) = c(proc$more$parnames, proc$more$betanames)
    if(!is.null(lik$report)){ print(f) }
    f = sgn*f
    df = sgn*df
    return(list(f = f, df = df, coefs = ncoefs))
}



SplineCoefsDC.DDE <- function(coefs, times, data, lik, proc, pars, beta, sgn = 1, basisvals, fdobj0){
    fdnames <- list(NULL, NULL, NULL)
    coefs2 = matrix(coefs, ncol(lik$bvals), length(coefs)/ncol(lik$bvals))
    fdnames[[2]] <- attr(coefs2, "dimnames")[[2]]
    fdobj.d <- list(coefs = coefs2, basis = basisvals, fdnames =fdnames)
    attr(fdobj.d, "class") <- "fd"
    delayProcObj <- delay.fit.sparse(fd0 = fdobj0, fd.d = fdobj.d, times = proc$more$qpts, tau = proc$more$more$tau, beta = beta, ndelay = proc$more$more$ndelay )
    delayLikObj <- delay.fit.sparse(fd0 = fdobj0, fd.d = fdobj.d, times = times,tau = lik$more$more$tau, beta= beta, ndelay = lik$more$more$ndelay)
    lik$more$more$y.d <- delayLikObj$y.d
    proc$more$more$y.d <- delayProcObj$y.d
    lik$more$more$bvals.d <- delayLikObj$bvals.d
    proc$more$more$bvals.d <- delayProcObj$bvals.d
    proc$more$more$bvals.d.list <- delayProcObj$bvals.d.list
    proc$more$more$y.d.list <- delayProcObj$y.d.list
    devals = as.matrix(lik$bvals %*% coefs2)
    colnames(devals) = proc$more$names
    g = as.matrix(t(lik$bvals) %*% lik$dfdx(data, times, devals,
        pars, lik$more)) + proc$dfdc(coefs2, proc$bvals, pars,
        proc$more)
    g = as.vector(g)
    return(sgn * g)
}

SplineCoefsDC2.DDE <- function(coefs, times, data, lik, proc, pars, beta, sgn = 1, basisvals, fdobj0){
    fdnames <- list(NULL, NULL, NULL)
    coefs2 = matrix(coefs, ncol(lik$bvals), length(coefs)/ncol(lik$bvals))
    fdnames[[2]] <- attr(coefs2, "dimnames")[[2]]
    fdobj.d <- list(coefs = coefs2, basis = basisvals, fdnames =fdnames)
    attr(fdobj.d, "class") <- "fd"
    delayProcObj <- delay.fit.sparse(fd0 = fdobj0, fd.d = fdobj.d, times = proc$more$qpts, tau = proc$more$more$tau, beta = beta, ndelay = proc$more$more$ndelay )
    delayLikObj <- delay.fit.sparse(fd0 = fdobj0, fd.d = fdobj.d, times = times,tau = lik$more$more$tau, beta= beta, ndelay = lik$more$more$ndelay)
    lik$more$more$y.d <- delayLikObj$y.d
    proc$more$more$y.d <- delayProcObj$y.d
    lik$more$more$bvals.d <- delayLikObj$bvals.d
    proc$more$more$bvals.d <- delayProcObj$bvals.d
    proc$more$more$bvals.d.list <- delayProcObj$bvals.d.list
    proc$more$more$y.d.list <- delayProcObj$y.d.list
    result = as.matrix(SplineCoefsDC2sparse(coefs, times, data, lik, proc, pars, sgn))
    return(result)
}

SplineCoefsErr.DDE <- function(coefs, times, data, lik, proc, pars, beta, sgn = 1, basisvals, fdobj0){
    fdnames <- list(NULL, NULL, NULL)
    coefs2 = matrix(coefs, ncol(lik$bvals), length(coefs)/ncol(lik$bvals))
    fdnames[[2]] <- attr(coefs2, "dimnames")[[2]]
    fdobj.d <- list(coefs = coefs2, basis = basisvals, fdnames =fdnames)
    attr(fdobj.d, "class") <- "fd"
    delayProcObj <- delay.fit.sparse(fd0 = fdobj0, fd.d = fdobj.d, times = proc$more$qpts, tau = proc$more$more$tau, beta = beta, ndelay = proc$more$more$ndelay )
    delayLikObj <- delay.fit.sparse(fd0 = fdobj0, fd.d = fdobj.d, times = times,tau = lik$more$more$tau, beta= beta, ndelay = lik$more$more$ndelay, lik = TRUE)
    lik$more$more$y.d <- delayLikObj$y.d
    proc$more$more$y.d <- delayProcObj$y.d
    ## lik$more$more$bvals.d <- delayLikObj$bvals.d
    proc$more$more$bvals.d <- delayProcObj$bvals.d
    proc$more$more$bvals.d.list <- delayProcObj$bvals.d.list
    proc$more$more$y.d.list <- delayProcObj$y.d.list
    devals = as.matrix(lik$bvals %*% coefs2)
    colnames(devals) = proc$more$names
    f = sum(lik$fn(data, times, devals, pars, lik$more)) + proc$fn(coefs2,
        proc$bvals, pars, proc$more)
    if (!is.null(proc$report)) {
        print(f)
    }
    return(sgn * f)
}

SplineCoefsDCDP.sparse <- function (coefs, times, data, lik, proc, pars, sgn = 1)
{
    coefs2 = matrix(coefs, ncol(lik$bvals), length(coefs)/ncol(lik$bvals))
    ## Till now, H has all 0's
    H <- proc$d2fdcdp(coefs2, proc$bvals, pars, proc$more)
    return(as.matrix(sgn * H))
}

LS.sparse <- function(fn, data, times, basisvals = NULL,
    lambda, fd.obj = NULL, more = NULL, weights = NULL, quadrature = NULL,
    in.meth = "nlminb", out.meth = "nls", control.in = list(),
    control.out = list(), eps = 1e-06, active = NULL, posproc = FALSE,
    poslik = FALSE, discrete = FALSE, names = NULL, sparse = FALSE,
    likfn = make.id(), likmore = NULL, delay = NULL, tauMax = NULL,
    basisvals0 = NULL, coefs0 = NULL, nbeta, ndelay, tau, nnls.res)
{
    ##################################################
    ## Added to orginal Profile.LS
    ## if(is.null(tauMax)) tauMax <- min(times) + 1/3 * (range(times)[2]- range(times)[1])
    ## Prepare to calculate the delay:
    ## data.d <- data[times>=tauMax, ,drop = FALSE]
    ## times.d <- knots.d <- times[times >= tauMax]
    ## norder <- basisvals$nbasis - length(basisvals$params)
    ## nbasis.d <- length(knots.d) + norder - 2
    ## range.d <- range(knots.d)
    ## basisvals.d and basisvals.0 are supplied rather than created
    ## basisvals.d <- create.bspline.basis(range=range.d, nbasis=nbasis.d, norder=norder, breaks=knots.d)
    betanames <- c()
    for(i in 1:length(nbeta)){
        for(j in 1:nbeta[i]){
            betanames <- c(betanames,paste("beta",i,".",j, sep = ""))
        }
    }
    pars <- nnls.res$pars
    coefs <- nnls.res$coefs
    beta <- nnls.res$beta
    if (is.null(active)) {
        active = 1:length(pars)
    }
    apars <- pars[active]
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
    delayProcObj <- delay.fit.sparse(fd0 = fdobj0, fd.d = fdobj.d, times = proc$more$qpts, tau = tau, beta= beta, ndelay = ndelay )
    delayLikObj <- delay.fit.sparse(fd0 = fdobj0, fd.d = fdobj.d, times = times,tau = tau, beta= beta, ndelay = ndelay)
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
    aparamnames = names(apars)
    if (is.null(control.out$maxIter)) {
        control.out$maxIter = 100
    }
    if (is.null(control.out$tol)){
        control.out$tol = 1e-08
    }
    res <- nnls.res
    ncoefs <- res$coefs
    if(is.null(control.out$pars.c))
        control.out$pars.c <- 100
    if(control.out$lambda.sparse == -1){
        Zdf <- res$Zdf
        Xdf <- res$Xdf
        y <- res$y
        lambda0 = max(abs(as.vector(t(y) %*% Zdf)))
        lambda = exp(seq(log(lambda0), log(lambda0 * 0.001), len = 20))
        pars.pen <- beta.pen <- coefs.pen <- list()
        bic <- f <- rep(NA, length(lambda))
        for(i in 1:length(lambda)){
            lambda.sparse <- lambda[i]
            res.sparse <- penalized(response = y, penalized = Zdf, unpenalized = Xdf, lambda1 = lambda[i], positive = TRUE, trace = FALSE)
            pars.pen[[i]] <- res.sparse@unpenalized
            beta.pen[[i]] <- res.sparse@penalized
            Ires <- inneropt.DDE(data, times, par = pars.pen[[i]], beta = beta.pen[[i]],  ncoefs, lik, proc, in.meth, control.in, basisvals = basisvals, fdobj0 = fdobj0)
            devals <- as.matrix(lik$bvals%*%Ires$coefs)
            f <- as.vector(as.matrix(data - lik$more$fn(times, devals, pars, lik$more$more)))
            coefs.pen[[i]] <- Ires$coefs
            sd.pen <- sd(f)
            ll.pen <- - sum(f^2) / (sd.pen^2) / 2 - length(f) * log(sd.pen)
            bic[i] <- -2 * ll.pen + (sum(beta.pen[[i]] > .Machine$double.eps) + length(pars.pen[[i]])) * log(length(data))
        }
        i.select <- which(bic == min(bic))
        sel.res <- list(pars.pen = pars.pen[[i.select]], beta.pen = beta.pen[[i.select]], bic = bic[i.select], coefs.pen = coefs.pen[[i.select]], lambda = lambda[i.select])
    }

    if(control.out$lambda.sparse == -2){
        ## Positive addaptive lasso using lars: Total eliminatioin
        ## NOT WORKING!!
        w.pars <- abs(res$pars) * control.out$pars.c
        w.beta <- abs(res$beta)                      # weights for adaptive lasso
        beta.ind <- which(w.beta > 0)
        w.beta <- w.beta[beta.ind]
        w <- c(w.pars, w.beta)
        x <- cbind(res$Xdf, res$Zdf[, beta.ind])
        n <- nrow(x)
        one <- rep(1, n)
        meanx <- drop(one %*% x)/n
        xc <- scale(x, meanx, FALSE)         # first subtracts mean
        normx <- sqrt(drop(one %*% (xc^2)))
        names(normx) <- NULL
        xs <- scale(xc, FALSE, normx)        # now rescales with norm (not sd)
        xs <- scale(xs, center=FALSE,scale=1/w)  # xs times the weights
        object <- lars.pos(x = xs, y = res$y, positive = TRUE, normalize = FALSE)
        object$beta <- object$beta[which(rowSums(object$beta[,1:length(res$pars), drop=FALSE] > 0) > 0 ),, drop = FALSE]
        object$beta <- sweep(object$beta, 1, w, "*")
        bic <- rep(0, dim(object$beta)[2])
        coefs.pen <- list()
        for(i in 1:dim(object$beta)[1]){
            pars.al1 <- object$beta[1:length(pars),i]
            names(pars.al1) <- names(pars)
            beta.al1 <- rep(0, length(beta))
            beta.al1[beta.ind] <- object$beta[(1+length(pars)) : dim(object$beta)[1],i]
            Ires <- inneropt.DDE(data, times, par = pars.al1, beta = beta.al1,  ncoefs, lik, proc, in.meth, control.in, basisvals = basisvals, fdobj0 = fdobj0)
            devals <- as.matrix(lik$bvals%*%Ires$coefs)
            f <- as.vector(as.matrix(data - lik$more$fn(times, devals, pars, lik$more$more)))
            coefs.pen[[i]] <- Ires$coefs
            sd.pen <- sd(f)
            ll.pen <- - sum(f^2) / (sd.pen^2) / 2 - length(f) * log(sd.pen)
            bic[i] <- -2 * ll.pen + (sum(beta.al1 > .Machine$double.eps) + length(pars.al1)) * log(length(data))
        }
        i.select <- which(bic == min(bic))
        pars.al1 <- object$beta[1:length(pars),i.select]
        names(pars.al1) <- names(pars)
        beta.al1 <- rep(0, length(beta))
        beta.al1[beta.ind] <- object$beta[(1+length(pars)) : dim(object$beta)[1],i.select]
        sel.res <- list(pars = pars.al1, beta = beta.al1, bic = bic[i.select], coefs = coefs.pen[[i.select]], lambda = lambda[i.select])
    }
    ## Positive addaptive lasso using lars: Partial eliminatioin:
    if(control.out$lambda.sparse == -3){
        y <- res$y - res$Xdf %*% res$pars
        w.beta <- abs(res$beta)                      # weights for adaptive lasso
        w.beta[w.beta == 0] <- min(w.beta[w.beta > 0]) / 2
        x <- res$Zdf
        n <- nrow(x)
        one <- rep(1, n)
        meanx <- drop(one %*% x)/n
        xc <- scale(x, meanx, FALSE)         # first subtracts mean
        normx <- sqrt(drop(one %*% (xc^2)))
        names(normx) <- NULL
        xs <- scale(xc, FALSE, normx)        # now rescales with norm (not sd)
        xs <- scale(xs, center=FALSE, scale=1/w.beta)  # xs times the weights
        object <- lars.pos(xs, y, type="lasso",normalize=FALSE, positive = TRUE)
        object$beta <- sweep(object$beta, 2, w.beta / normx, "*")
        bic <- rep(0, dim(object$beta)[1])
        coefs.pen <- list()
        for(i in 1:dim(object$beta)[1]){
            beta.al <- object$beta[i,]
            Ires <- inneropt.DDE(data, times, par = res$pars, beta = beta.al,  ncoefs, lik, proc, in.meth, control.in, basisvals = basisvals, fdobj0 = fdobj0)
            devals <- as.matrix(lik$bvals%*%Ires$coefs)
            f <- as.vector(as.matrix(data - lik$more$fn(times, devals, res$pars, lik$more$more)))
            coefs.pen[[i]] <- Ires$coefs
            sd.pen <- sd(f)
            ll.pen <- - sum(f^2) / (sd.pen^2) / 2 - length(f) * log(sd.pen)
            bic[i] <- -2 * ll.pen + (sum(beta.al > .Machine$double.eps) + length(res$pars)) * log(length(data))
        }
        i.select <- which(bic == min(bic))
        beta.al <- object$beta[i.select,]
        sel.res <- list(pars = res$pars, beta = beta.al, bic = bic[i.select], coefs = coefs.pen[[i.select]])
    }
    ## Positive Lars:
    if(control.out$lambda.sparse == -4){
        Zdf <- res$Zdf
        y <- res$y - res$Xdf %*% res$pars
        object <- lars.pos(Zdf,y, positive = TRUE)
        bic <- rep(0, dim(object$beta)[1])
        coefs.pen <- list()
        for(i in 1:dim(object$beta)[1]){
            beta.pl <- object$beta[i,]
            Ires <- inneropt.DDE(data, times, par = res$pars, beta = beta.pl,  ncoefs, lik, proc, in.meth, control.in, basisvals = basisvals, fdobj0 = fdobj0)
            devals <- as.matrix(lik$bvals%*%Ires$coefs)
            f <- as.vector(as.matrix(data - lik$more$fn(times, devals, pars, lik$more$more)))
            coefs.pen[[i]] <- Ires$coefs
            sd.pen <- sd(f)
            ll.pen <- - sum(f^2) / (sd.pen^2) / 2 - length(f) * log(sd.pen)
            bic[i] <- -2 * ll.pen + (sum(beta.pl > .Machine$double.eps) +length(res$pars)) * log(length(data))
        }
        i.select <- which(bic == min(bic))
        beta.al <- object$beta[i.select ,]
        sel.res <- list(pars = res$pars, beta = beta.al, bic = bic[i.select], coefs = coefs.pen[[i.select]])
    }
    if(control.out$lambda.sparse < 0){
        return(list( data = data,res = res, select = sel.res))
    }else{
        return(list( data = data,res = res))
    }

}
