library(CollocInfer)
library(deSolve)

blowfly.gen <- function(t, y, parms){
    if(t<0)
        lag <- 1000
    else
        lag <- lagvalue(t - parms["tau"])
    dy <- parms["c"] * lag * exp(-lag / parms["N0"]) - parms["a"] * y
    list(dy, dy = dy)
}

blowfly.pars <- c(150 / 8, 8 / 8 , 1000, 8)
names(blowfly.pars) <- c("c", "a", "N0", "tau")
times <- seq(-10, 200, by = 0.5)
yout <- dede(y = 1000, times = times, func = blowfly.gen, parms = blowfly.pars, atol = 1e-7)
## plot(yout[55:405,], which = 1, type = "l", lwd = 2, main = "Nicholson's Blowflies Model")

set.seed(42)
data.res <- list()
for(i in 1:500){
    data.res[[i]] <- list()
    blowfly.data <- yout[,2][55:405]   + rnorm(351, sd = 250)
    initPars <- blowfly.pars[1:3] + c(runif(1, -4, 4), runif(1, -0.075, 0.075), runif(1, -200, 200))
    names(initPars) <- names(blowfly.pars)[1:3]
    data.res[[i]]$blowfly.data <- blowfly.data
    data.res[[i]]$initPars <- initPars
    initBeta <- rep(1/10, 10) + runif(10, -0.02, 0.02)
    data.res[[i]]$initBeta <- initBeta / sum(initBeta)
}

curseed <- get(".Random.seed", .GlobalEnv)
save(data.res, blowfly.pars, curseed, file = "data-blowfly-250.RData")
