Profile.LS.tv.delay <- function(fn, data, times, pars, beta, kappa, coefs = NULL, basisvals = NULL,
    lambda, fd.obj = NULL, more = NULL, weights = NULL, quadrature = NULL,
    in.meth = "nlminb", out.meth = "nls", control.in = list(),
    control.out = list(), eps = 1e-06, active = NULL, posproc = FALSE,
    poslik = FALSE, discrete = FALSE, names = NULL, sparse = FALSE,
    likfn = make.id(), likmore = NULL, delay = NULL, tauMax = NULL,
    basisvals0 = NULL, coefs0 = NULL, nbeta, ndelay, tau)
{
    if (is.null(active)) {
        active = 1:length(c(pars, kappa))
    }

    ## Create y.d
    fdnames <- list(NULL, NULL, NULL)
    fdnames[[2]] <- attr(coefs, "dimnames")[[2]]
    fdobj0 <- list(coefs = coefs0, basis = basisvals0, fdnames =fdnames)
    fdobj.d <- list(coefs = coefs, basis = basisvals, fdnames =fdnames)
    attr(fdobj0, "class") <- "fd"
    attr(fdobj.d, "class") <- "fd"

    profile.obj = LS.setup(pars = c(pars, kappa), coefs = coefs, fn = fn,
    basisvals, lambda = lambda, fd.obj, more, data, weights,
        times, quadrature, eps = 1e-06, posproc, poslik, discrete,
        names, sparse, likfn = make.id(), likmore = NULL)
    dims = dim(data)
    lik = profile.obj$lik
    proc = profile.obj$proc
    proc$more$more$nKappa <- length(kappa)
    coefs = profile.obj$coefs
    data = profile.obj$data
    times = profile.obj$times

    ## Create names for delay parameters beta
    betanames <- c()
    for(i in 1:length(nbeta)){
        for(j in 1:nbeta[i]){
            betanames <- c(betanames,paste("beta",i,".",j, sep = ""))
        }
    }
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


    Ires <- inneropt.DDE(data, times, c(pars, kappa), beta, coefs, lik, proc, in.meth, control.in, basisvals = basisvals, fdobj0 = fdobj0)
    ## Ires <- IresTmp
    ncoefs <- Ires$coefs
    apars = c(pars, kappa)[active]
    aparamnames = names(apars)
    if (is.null(control.out$maxIter)) {
        control.out$maxIter = 100
    }
    if (is.null(control.out$tol)){
        control.out$tol = 1e-08
    }
     if (is.null(control.out$echo)){
         control.out$echo = TRUE
    }
    res <- nls.tv.delay(pars = pars, beta = beta, kappa = kappa, active = active, basisvals = basisvals, fdobj0 = fdobj0, times = times, data = data, coefs = ncoefs, lik = lik, proc = proc, control.out = control.out, control.in = control.in, in.meth = in.meth)
    ncoefs <- res$coefs
    return(list( data = data,res = res, ncoefs = ncoefs))
}


nls.tv.delay <- function(pars, beta, kappa, active, basisvals, fdobj0, times, data, coefs, lik, proc, start, X.index, control.out, control.in, in.meth){
    if(control.out$method == "twoStage"){
        delta <- rep(1, length(data))
    }
    pars.names <- names(pars)
    kappa.names <- names(kappa)
    f.conv <- pars.kappa.beta <- c()
    maxStep <- 10
    lambda.sparse <- control.out$lambda.sparse
    for(i in 1:control.out$maxIter){
        for(j in 1:maxStep){
            linObj <- ProfileSSE.AllPar.sparse(pars = c(pars,kappa), beta = beta, times = times, data = data, coefs = coefs, lik = lik, proc = proc, in.meth = in.meth,control.in = control.in, basisvals = basisvals, fdobj0 = fdobj0)
            f.new <- linObj$f
            f.new <- sum(f.new^2)
            if(control.out$echo){
                print(x = c(paste("Iter:", i, f.new)))
                cat(pars, kappa, beta, "\n")
            }
            if(i == 1){
                break
            }else{
                if(f.conv[i - 1] - f.new > 0 & f.conv[i - 1] - f.new < control.out$tol){
                    return(list(pars=pars.old, kappa = kappa.old, coefs = coefs, f = f.new, y = y, Xdf = Xdf, Zdf = Zdf, conv = list(f = f.conv, pars.kappa.beta=pars.kappa.beta, conv.message = "Converged.")))
                }
                if(f.conv[i - 1] - f.new > 0){
                    break
                }
                if(f.conv[i - 1] - f.new < 0 & j == maxStep){
                    return(list(pars=pars.old, kappa = kappa.old, coefs = coefs, f = f.new,  y = y, Xdf = Xdf, Zdf = Zdf, conv = list(f = f.conv, pars.kappa.beta=pars.kappa.beta, conv.message = "Non-dereasing objective.")))
                }
                pars <- 0.5*(pars - pars.old) + pars.old
                kappa <- 0.5*(kappa - kappa.old) + kappa.old
                beta <- 0.5*(beta - beta.old) + beta.old
            }
        }
        pars.old <- pars
        kappa.old <- kappa
        beta.old <- beta
        f.conv <- c(f.conv, f.new)
        pars.kappa.beta <- rbind(pars.kappa.beta, c(pars, kappa, beta))
        Xdf <- - linObj$df[, 1:length(pars), drop = FALSE]
        Zdf <- - linObj$df[, (length(pars) +1): dim(linObj$df)[2]]
        y <- - linObj$df %*% c(pars, kappa, beta) + linObj$f
        coefs <- linObj$coefs
        if(control.out$method == "nnls"){
            E <- t(c(rep(0, length(pars)) , rep(0, length(kappa)), rep(1, length(beta))))
            F <- 1
            G <- diag(length(c(pars, kappa, beta)))
            H <- rep(0, length(c(pars, kappa, beta)))
            res <- lsei(A= cbind(Xdf, Zdf), B = y, E = E, F=F, G = G, H = H)
            ## res <- nnls(A = cbind(Xdf, Zdf), b= y)
            pars <- res$X[1:length(pars)]
            kappa <- res$X[(length(pars) + 1) : (length(pars) + length(kappa))]
            beta <- res$X[(length(pars) + length(kappa) + 1) : (length(pars) + length(kappa) + length(beta))]
        }
        names(pars) <- pars.names
        names(kappa) <- kappa.names
    }
    return(list(pars=pars.old, kappa = kappa.old, beta = beta.old, coefs = coefs, f = f.new, y = y, Xdf = Xdf, Zdf = Zdf, conv = list(f = f.conv, pars.kappa.beta=pars.kappa.beta,  conv.message = "Maximum iterations reached.")))
}


sparse.tv.delay <- function(fn, data, times, pars, beta, kappa, coefs = NULL, basisvals = NULL,
    lambda, fd.obj = NULL, more = NULL, weights = NULL, quadrature = NULL,
    in.meth = "nlminb", out.meth = "nls", control.in = list(),
    control.out = list(), eps = 1e-06, active = NULL, posproc = FALSE,
    poslik = FALSE, discrete = FALSE, names = NULL, sparse = FALSE,
    likfn = make.id(), likmore = NULL,delay, basisvals0 = NULL, coefs0 = NULL,  nbeta, ndelay, tau, nnls.res)
{
    betanames <- c()
    for(i in 1:length(nbeta)){
        for(j in 1:nbeta[i]){
            betanames <- c(betanames,paste("beta",i,".",j, sep = ""))
        }
    }

    ## pars <- nnls.res$pars
    ## beta <- nnls.res$beta
    ## kappa <- nnls.res$kappa
    coefs <- nnls.res$coefs

    if (is.null(active)) {
        active = 1:length(c(pars,beta, kappa))
    }
    apars <- pars[active]

    ## Create y.d
    fdnames <- list(NULL, NULL, NULL)
    fdnames[[2]] <- attr(coefs, "dimnames")[[2]]
    fdobj0 <- list(coefs = coefs0, basis = basisvals0, fdnames =fdnames)
    fdobj.d <- list(coefs = coefs, basis = basisvals, fdnames =fdnames)
    attr(fdobj0, "class") <- "fd"
    attr(fdobj.d, "class") <- "fd"

    profile.obj = LS.setup(pars = c(pars, kappa), coefs = coefs, fn = fn,
    basisvals, lambda = lambda, fd.obj, more, data, weights,
        times, quadrature, eps = 1e-06, posproc, poslik, discrete,
        names, sparse, likfn = make.id(), likmore = NULL)
    dims = dim(data)
    lik = profile.obj$lik
    proc = profile.obj$proc
    proc$more$more$nKappa <- length(kappa)
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
        if(is.null(control.out$maxInnerIter)) maxInnerIter <- 10
        if(is.null(control.out$lambda.len)) nlambda <- 3

        Zdf1 <- res$Zdf[, 1: length(kappa), drop = FALSE]
        Zdf2 <- res$Zdf[,(1+length(kappa)):(length(kappa)+length(beta)),drop = FALSE]
        Xdf <- res$Xdf

        y <- res$y
        Z1 <- rowSums(Zdf1)
        Z2 <- rowSums(Zdf2)

        beta0 <- sum((y- Xdf %*% pars) * Z1) / sum(Z1^2)
        lambda10 <- max(abs(as.vector(t(y -  Xdf %*% pars - Z1 * beta0) %*% sweep(Zdf1, 1, mean(Zdf1)))))
        lambda1 = exp(seq(log(lambda10), log(lambda10 * 0.01), len = nlambda))
        lambda20 <- max(abs(as.vector(t(y) %*% Zdf2)))
        lambda2 = exp(seq(log(lambda20), log(lambda20 * 0.01), len = nlambda))

        pars.pen <- kappa.pen <- beta.pen <- coefs.pen <- list()
        bic <- rep(NA, length(lambda1) * length(lambda2))
        for(i in 1:nlambda){
            pars.pen[[i]] <- kappa.pen[[i]] <- beta.pen[[i]] <- coefs.pen[[i]] <- list()
            for(j in 1:nlambda){
                ij <- (i - 1) * nlambda + j
                lambda.i1 <- lambda1[i]
                lambda.j1 <- lambda2[j]
                y1 <- y - Zdf2 %*% beta
                y2 <- y - Zdf1 %*% kappa
                kappa.pen[[ij]] <- kappa
                pars.pen[[ij]] <- pars
                beta.pen[[ij]] <- beta
                for(k in 1:maxInnerIter){
                    res.sparse1 <- penalized(response = y1, penalized = Zdf1, unpenalized = Xdf, lambda2 = lambda1[i], fusedl = TRUE, positive = TRUE, trace = FALSE)
                    res.sparse2 <- penalized(response = y2, penalized = Zdf2, unpenalized = Xdf, lambda1 = lambda2[j], positive = TRUE, trace = FALSE)
                    if(sum((kappa.pen[[ij]] - res.sparse1@penalized)^2, (pars.pen[[ij]] - (res.sparse1@unpenalized + res.sparse2@unpenalized)/2)^2, (res.sparse1@unpenalized  - res.sparse2@unpenalized)^2, (beta.pen[[ij]] - res.sparse2@penalized)^2) < eps){
                        kappa.pen[[ij]] <- res.sparse1@penalized
                        pars.pen[[ij]] <- (res.sparse2@unpenalized + res.sparse1@unpenalized) / 2
                        beta.pen[[ij]] <- res.sparse2@penalized
                        break
                    }
                    kappa.pen[[ij]] <- res.sparse1@penalized
                    pars.pen[[ij]] <- (res.sparse2@unpenalized + res.sparse1@unpenalized) / 2
                    beta.pen[[ij]] <- res.sparse2@penalized
                    y1 <- y - Zdf2 %*% beta.pen[[ij]]
                    y2 <- y - Zdf1 %*% kappa.pen[[ij]]
                }
                Ires <- inneropt.DDE(data, times, c(pars.pen[[ij]],kappa.pen[[ij]]), beta.pen[[ij]], coefs, lik, proc, in.meth, control.in, basisvals = basisvals, fdobj0 = fdobj0)
                devals <- as.matrix(lik$bvals%*%Ires$coefs)
                f <- as.vector(as.matrix(data - lik$more$fn(times, devals, pars, lik$more$more)))
                coefs.pen[[ij]] <- Ires$coefs
                sd.pen <- sd(f)
                ll.pen <- - sum(f^2) / (sd.pen^2) / 2 - length(f) * log(sd.pen)
                bic[(i-1) * length(lambda1) + j] <- -2 * ll.pen + (length(unique(kappa.pen[[ij]])) + length(pars.pen[[ij]]) + sum(beta.pen[[ij]] != 0)) * log(length(data))
            }
        }
        ij.select <- which(bic == min(bic))
        sel.res <- list(pars.pen = pars.pen[[ij.select]], kappa.pen = kappa.pen[[ij.select]], bic = bic[ij.select], coefs.pen = coefs.pen[[ij.select]], lambda = c(lambda1[ceiling(ij.select / nlambda)], lambda2[ifelse(ij.select %% nlambda == 0, nlambda, ij.select %% nlambda)]))
    }
    return(list(select = sel.res))
}
