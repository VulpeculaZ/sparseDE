## Fit nls model for blowfly model
source("./R/blowflies.R")
library(deSolve)
library(FME)

blowfly_nls <- function(pars, init){
    times <- seq(0,175, 0.5)
    yout <- dede(y = 1000, times = times, func = blowfly.gen, parms = pars, atol = 1e-7)
    as.data.frame(yout)
}

blowfly.pars <- c(150 / 8, 8 / 8 )
names(blowfly.pars) <- c("c", "a", "N0")

blowfly.fn <- function(t, y, parms){
    if(t<0)
        lag <- 1000
    else
        lag <- lagvalue(t - 8)
    dy <- parms["c"] * lag * exp(-lag / parms["N0"]) - parms["a"] * y
    list(dy, dy = dy)
}

blowfly.gen <- function(pars){
    names(pars) <- c("c", "a", "N0")
    times <- seq(-10, 200, by = 0.5)
    yout <- dede(y = 1000, times = times, func = blowfly.fn, parms = pars, atol = 1e-7)
    as.data.frame(yout)
    colnames(yout) <- c("time","y","dy.c")
    return(yout)
}

blowfly.cost <- function(pars){
    out <- blowfly.gen(pars)
    cost <- modCost(model = out, obs = blowfly.data, err = "sd")
    return(cost)
}


Blowfly.pars <- c(150 / 8, 8 / 8 , 1000)
blowfly.pars <- blowfly.pars + runif(3, -0.1, 0.1) * blowfly.pars
names(blowfly.pars) <- c("c", "a", "N0")
times <- seq(-10, 200, by = 0.5)
yout <- dede(y = 1000, times = times, func = blowfly.fn, parms = blowfly.pars, atol = 1e-7)
## plot(yout[55:405,], which = 1, type = "l", lwd = 2, main = "Nicholson's Blowflies Model")

set.seed(42)
blowfly.data <- yout
blowfly.data[,2] <- blowfly.data[,2] + rnorm(421, sd = 250)
blowfly.data <- cbind(blowfly.data, rep(250, 421))
blowfly.data <- cbind(time = yout[,1],
             y = yout[,2] + rnorm(n = 421, sd = 250),
             sd = 250)
blowfly.cost(blowfly.pars)$model
blowfly.fit <- modFit(f = blowfly.cost, p = blowfly.pars)


HIV_R <- function (pars, V_0 = 50000, dV_0 = -200750, T_0 = 100) {
    derivs <- function(time, y, pars){
        with (as.list(c(pars, y)), {
            dT <- lam - rho * T - bet * T * V
            dI <- bet * T * V - delt * I
            dV <- n * delt * I - c * V - bet * T * V
            return(list(c(dT, dI, dV), logV = log(V)))
        })
    }
    I_0   <- with(as.list(pars), (dV_0 + c * V_0) / (n * delt))
    y <-c(T=T_0,I=I_0,V=V_0)
    times <- c(seq(0, 0.8, 0.1), seq(2, 60, 2))
    out <- ode(y = y, parms = pars, times = times, func = derivs)
    as.data.frame(out)
}

HIVcost <- function(pars){
    out <- HIV(pars)
    cost <- modCost(model = out, obs = DataLogV, err = "sd")
    return(modCost(model = out, obs = DataT, err = "sd", cost = cost))
}



