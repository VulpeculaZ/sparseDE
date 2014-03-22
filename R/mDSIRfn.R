transm <- function(t){
    month <- t %% 1
    r <- rep(0.5, length(t))
    r[month  > 6/12 & month < 9 /12] <- 0
    return(r)
}

mDSIRfn <- list()

mDSIRfn$fn <- function (t, y, p, more)
{
    r = y
    yi.d <- more$y.d[,1]
    b <- more$b
    r[, "S"] =  - (p["sig"] + transm(t)) * yi.d * y[, "S"] + b ## p["alpha"]
    r[, "I"] =  (p["sig"] + transm(t))* yi.d * y[, "S"] - p["gamma"] * y[, "I"]
    return(r)
}

mDSIRfn$dfdx <- function (t, y, p, more)
{
    r = array(0, c(length(t), ncol(y), ncol(y)))
    dimnames(r) = list(NULL, colnames(y), colnames(y))
    yi.d <- more$y.d[,1]
    r[, "S", "S"] = - (p["sig"] + transm(t)) * yi.d
    r[, "I", "S"] =  (p["sig"] + transm(t)) * yi.d
    r[, "I", "I"] = -p["gamma"]
    return(r)
}

mDSIRfn$dfdx.d <- function (t, y, p, more)
{
    ## yi.d <- more$y.d[,2]
    r = array(0, c(length(t), ncol(y), ncol(y)))
    dimnames(r) = list(NULL, colnames(y), colnames(y))
    r[, "S", "I"] = -  (p["sig"] + transm(t)) * y[,"S"]
    r[, "I", "I"] =  (p["sig"] + transm(t)) * y[,"S"]
    return(r)
}

mDSIRfn$dfdp <- function (t, y, p, more)
{
    b <- more$b
    yi.d <- more$y.d[,1]
    r = array(0, c(length(t), ncol(y), length(p)))
    dimnames(r) = list(NULL, colnames(y), names(p))
    r[ , "S", "sig"] = - yi.d * y[, "S"]
    ## r[ , "S", "alpha"] = b
    r[, "I", "gamma"] = - y[, "I"]
    r[ , "I", "sig"] = yi.d * y[, "S"]
    return(r)
}

mDSIRfn$d2fdx2 <- function (t, y, p, more)
{
    r = array(0, c(length(t), ncol(y), ncol(y), ncol(y)))
    dimnames(r) = list(NULL, colnames(y), colnames(y), colnames(y))
    return(r)
}

mDSIRfn$d2fdxdp <- function (t, y, p, more)
{
    yi.d <- more$y.d[,1]
    r = array(0, c(length(t), ncol(y), ncol(y), length(p)))
    dimnames(r) = list(NULL, colnames(y), colnames(y), names(p))
    r[, "I", "I", "gamma"] = -1
    r[, "S", "S", "sig"] = -yi.d
    r[, "I", "S", "sig"] = yi.d
    return(r)
}


mDSIRfn$d2fdx.ddp <- function (t, y, p, more)
{
    yi.d <- more$y.d[,1]
    r = array(0, c(length(t), ncol(y), ncol(y), length(p)))
    dimnames(r) = list(NULL, colnames(y), colnames(y), names(p))
    r[, "S", "I", "sig"] = -y[,"S"]
    r[, "I", "I", "sig"] =  y[,"S"]
    return(r)
}

mDSIRfn$d2fdxdx.d <- function (t, y, p, more)
{
    yi.d <- more$y.d[,1]
    r = array(0, c(length(t), ncol(y), ncol(y), ncol(y)))
    dimnames(r) = list(NULL, colnames(y), colnames(y), colnames(y))
    r[, "S", "S", "I"] = -(p["sig"] + transm(t))
    r[, "I", "S", "I"] =  (p["sig"] + transm(t))
    return(r)
}


mDSIRfn$d2fdx.d2 <- function (t, y, p, more)
{
    yi.d <- more$y.d[,1]
    r = array(0, c(length(t), ncol(y), ncol(y), ncol(y)))
    dimnames(r) = list(NULL, colnames(y), colnames(y), colnames(y))
    return(r)
}
