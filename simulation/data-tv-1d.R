library(deSolve)


## Function to simulate DSIR
dtvSIR.gen <- function(t, y, parms, kappa){
    if(t < 0){
        lagI1 <- 2000
        ## lagI2 <- 400
    }
    else{
        lagI1 <- lagvalue(t - parms["tau1"], 2)
    }
    sint <- 8000 * (sin(t / parms["f"] / pi) / 2 + 2)
    dyS <- - (tvtrans(t, kappa)) * (lagI1) * y[1] +  sint
    dyI <- (tvtrans(t, kappa)) * (lagI1) * y[1] - parms["gamma"] * y[2]
    list(c(dyS, dyI))
}


yinit <- c(4000, 2000)
tvSIR.pars <- c(10, 1)
tau <- c(8/52)
dtvSIR.pars <- c(tvSIR.pars, tau)
names(dtvSIR.pars) <- c("gamma", "f", "tau1")
times <- seq(-dtvSIR.pars["tau1"], 5, by = 1/52)
kappa <- c(rep(0.005, 6), rep( 0.0025, 6))
yout <- dede(y = yinit, times = times, func = dtvSIR.gen, parms = dtvSIR.pars, atol = 1e-10, kappa = kappa)


set.seed(42)
data.res <- list()
for(i in 1:500){
    xout <- c()
    data.res[[i]] <- list()
    xout <- cbind(xout, yout[,2] + rnorm(length(yout[,2]), sd = 100))
    xout <- cbind(xout, yout[,3] + rnorm(length(yout[,2]), sd = 100))
    ## points(times, xout)
    xout0 <- xout[times >= 0,]
    data.res[[i]]$xout <- xout
    initPars <- 10 + runif(1, -1, 1)
    names(initPars) <- "gamma"
    data.res[[i]]$initPars <- initPars
    initBeta <- rep(1/6, 6) + runif(6, -0.02, 0.02)
    data.res[[i]]$initBeta <- initBeta
    initKappa <- kappa + runif(6, -0.001, 0.001)
    names(initKappa) <- paste("k", 1:length(initKappa), sep = "")
    data.res[[i]]$initKappa <- initKappa
}

curseed <- get(".Random.seed", .GlobalEnv)
save(data.res, dtvDSIR.pars, curseed, file = "data-tv-1d-sd100.RData")
