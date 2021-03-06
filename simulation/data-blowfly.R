library(CollocInfer)
library(deSolve)

blowfly.gen <- function(t, y, parms){
    lag <- rep(1000, 10)
    for(i in 1:10){
        if(t>0)
            lag[i] <- lagvalue(t - 8)
    }
    lag <- mean(lag)
    dy <- parms["c"] * lag * exp(-lag / parms["N0"]) - parms["a"] * y
    list(dy, dy = dy)
}


blowfly.pars <- c(150 / 8, 8 / 8 , 1000)
names(blowfly.pars) <- c("c", "a", "N0")
blowfly.pars <- c(blowfly.pars, seq(7.75, 8.25, length.out = 11))
times <- seq(-10, 200, by = 0.5)
yout <- dede(y = 1000, times = times, func = blowfly.gen, parms = blowfly.pars, atol = 1e-7)



set.seed(42)
data.res <- list()
for(i in 1:500){
    data.res[[i]] <- list()
    blowfly.data <- yout[,2][55:405]   + rnorm(351, sd = 250)
    ## pdf(file = "blowfly-sim.pdf", width = 7, height = 5)
    ## plot(x = times[55:405], y = yout[,2][55:405], type= "l", ylim = c(-500,8000), xlab = "t (day)", ylab = "y")
    ## points(x = times[55:405], y = blowfly.data)
    ## dev.off()
    initPars <- blowfly.pars[1:3] + c(runif(1, -4, 4), runif(1, -0.075, 0.075), runif(1, -200, 200))
    names(initPars) <- names(blowfly.pars)[1:3]
    data.res[[i]]$blowfly.data <- blowfly.data
    data.res[[i]]$initPars <- initPars
    initBeta <- rep(1/10, 10) + runif(10, -0.02, 0.02)
    data.res[[i]]$initBeta <- initBeta / sum(initBeta)
}

curseed <- get(".Random.seed", .GlobalEnv)
save(data.res, blowfly.pars, curseed, file = "data-blowfly-250-dist.RData")
