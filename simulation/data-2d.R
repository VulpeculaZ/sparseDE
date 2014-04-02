library(deSolve)


## Function to simulate DSIR
DSIR.gen <- function(t, y, parms){
    if(t < 0){
        lagI1 <- 0.05
        lagI2 <- 0.05
    }
    else{
        lagI1 <- lagvalue(t - parms["tau1"], 2)
        lagI2 <- lagvalue(t - parms["tau2"], 2)
    }
    sint <- sin(t / parms["f"]) + 1
    if(y[2] > 1)
        y[2] <- 1
    if(y[2] < 0)
        y[2] <- 0
    dyS <- -parms["beta"] * (lagI1 + lagI2) * y[1] + (1 - y[1]) * (parms["b"] * sint + 0.1)
    dyI <- parms["beta"] * (lagI1 + lagI2) * y[1] - (parms["gamma"] + (parms["b"] * sint + 0.1)) * y[2]
    list(c(dyS, dyI))
}


yinit <- c(0.9, 0.1)
times <- seq(0, 25, by = 0.1)
SIR.pars <- c(2, 0.25, 0.5, 1)
names(SIR.pars) <- c("beta", "b","gamma", "f")
## Simulation for DSIR
tau <- c(2,2+1/3)
DSIR.pars <- c(SIR.pars, tau)
names(DSIR.pars) <- c("beta", "b","gamma", "f", "tau1","tau2")
times <- seq(-DSIR.pars["tau2"], max(times), by = 0.1)
yout <- dede(y = yinit, times = times, func = DSIR.gen, parms = DSIR.pars, atol = 1e-10)
yout0 <- yout[times > 0, ]
times0 <- times[times >= 0]

set.seed(42)
data.res <- list()
for(i in 1:500){
    xout <- c()
    data.res[[i]] <- list()
    xout <- cbind(xout, yout[,2] + rnorm(length(yout[,2]), sd = 0.02))
    xout <- cbind(xout, yout[,3] + rnorm(length(yout[,2]), sd = 0.02))
    ## points(times, xout)
    xout0 <- xout[times >= 0,]
    data.res[[i]]$xout <- xout
    initPars <- c(0.5) + runif(1, -0.1, 0.1)
    names(initPars) <- c("gamma")
    data.res[[i]]$initPars <- initPars
    initBeta <- rep(0.1, 16) + runif(16, -0.02, 0.02)
    data.res[[i]]$initBeta <- initBeta
}
curseed <- get(".Random.seed", .GlobalEnv)
save(data.res, DSIR.pars, curseed, file = "data-2dadj-sd02.RData")