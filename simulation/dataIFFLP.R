library(deSolve)
tvI <- function(t,k){
    rt <- c(t[1], t[length(t)])
    stepLen <- (rt[2] -rt[1]) / length(k)
    r <- rep(NA, length(t))
    for(i in 1 : length(k)){
        r[t  >= rt[1] + (i-1) * stepLen & t < rt[1] + i * stepLen] <- k[i]
    }
    r[length(t)] <- k[length(k)]
    return(r)
}

Hill <- function(L, K){
    L / (L + K)
}


tvIFFLP.gen <- function(t, y, p){
    It <- 0.5
    if(t > 10) It <- 0.7
    dA <- It  * Hill(1 - y[1], p["KIA"]) - p["kFA"] * Hill(y[1], p["KFA"])
    dB <- y[1] * p["kAB"] * Hill(1 - y[1], p["KAB"]) - p["kFB"] * Hill(y[2], p["KFB"])
    dC <- y[1] * p["kAC"] * Hill(1 - y[3], p["KAC"]) - y[2] * p["kBC"] * Hill(y[3], p["KBC"])
    list(c(dA, dB, dC))
}

yinit <- c(0.2679492, 1.3559556, 0.1030671)
IFFLP.pars <- rep(1, 11)
names(IFFLP.pars) <- c("kFA", "kAB", "kFB", "kAC", "kBC", "KIA", "KFA", "KAB", "KFB", "KAC", "KBC")
IFFLP.pars["KFB"] <- 100
IFFLP.pars["KAB"] <- 0.001
IFFLP.pars["kAB"] <- 0.5
IFFLP.pars["kFB"] <- 10
kappa <- c(rep(0.5,2), rep(0.7, 6))
times <- seq(0, 40, by = 0.2)

yout <- dede(y = yinit, times = times, func = tvIFFLP.gen, parms = IFFLP.pars, atol = 1e-10)

## plot(yout)


set.seed(42)
data.res <- list()
for(i in 1:500){
    data.res[[i]] <- list()
    xout <- matrix(, nrow = length(yout[,1]), ncol = 3)
    xout[,1] <-  yout[,2] + rnorm(length(yout[,2]), sd = 0.01)
    xout[,2] <- yout[,3] + rnorm(length(yout[,3]), sd = 0.05)
    xout[,3] <- yout[,4] + rnorm(length(yout[,4]), sd = 0.003)
    ## points(times, xout)
    data.res[[i]]$xout <- xout
    initPars <- IFFLP.pars + IFFLP.pars * runif(length(IFFLP.pars), -0.3, 0.3)
    data.res[[i]]$initPars <- initPars
    initKappa <- kappa + runif(8, -0.1, 0.1)
    names(initKappa) <- paste("k", 1:length(initKappa), sep = "")
    data.res[[i]]$initKappa <- initKappa
    plot(yout[,1], yout[,4], type = "l")
    points(yout[,1], xout[,3])
}

curseed <- get(".Random.seed", .GlobalEnv)
save(data.res, IFFLP.pars, curseed, file = "data-IFFLP-01.RData")
