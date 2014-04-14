SIRfn <- list()

SIRfn$fn <- function (t, y, p, more)
{
    r = y
    r[, "S"] =  - p["beta"] * y[, "I"] * y[, "S"] + (1 - y[,"S"]) * (0.25 * (sin(t/2) + 1) + 0.1)
    r[, "I"] = p["beta"] *  y[, "I"] * y[, "S"] - (p["gamma"] + (0.25 * (sin(t/2) + 1) + 0.1)) * y[, "I"]
    return(r)
}

SIRfn$dfdp <- function (t, y, p, more)
{
    r = array(0, c(length(t), ncol(y), length(p)))
    dimnames(r) = list(NULL, colnames(y), names(p))
    r[, "S", "beta"] = - y[, "I"] * y[, "S"]
    r[, "I", "beta"] = y[, "I"] * y[, "S"]
    r[, "I", "gamma"] = - y[, "I"]
    return(r)
}

SIRfn$d2fdp2 <- function (t, y, p, more)
{
    r = array(0, c(length(t), ncol(y), length(p), length(p)))
    dimnames(r) = list(NULL, colnames(y), names(p), names(p))
    return(r)
}

SIRfn$dfdx <- function (t, y, p, more)
{
    r = array(0, c(length(t), ncol(y), ncol(y)))
    dimnames(r) = list(NULL, colnames(y), colnames(y))
    r[, "S", "S"] = - p["beta"] * y[, "I"] - (0.25 * (sin(t/2) + 1) + 0.1)
    r[, "S", "I"] = - p["beta"] * y[,"S"]
    r[, "I", "S"] = p["beta"] * y[, "I"]
    r[, "I", "I"] = p["beta"] * y[, "S"] - p["gamma"] - (0.25 * (sin(t/2) + 1) + 0.1)
    return(r)
}

SIRfn$d2fdx2 <- function (t, y, p, more)
{
    r = array(0, c(length(t), ncol(y), ncol(y), ncol(y)))
    dimnames(r) = list(NULL, colnames(y), colnames(y), colnames(y))
    r[, "S", "S", "I"] = - p["beta"]
    r[, "S", "I", "S"] = - p["beta"]
    r[, "I", "S", "I"] = p["beta"]
    r[, "I", "I", "S"] = p["beta"]
    return(r)
}

SIRfn$d2fdxdp <- function (t, y, p, more)
{
    r = array(0, c(length(t), ncol(y), ncol(y), length(p)))
    dimnames(r) = list(NULL, colnames(y), colnames(y), names(p))
    r[, "S", "S", "beta"] = - y[, "I"]
    r[, "S", "I", "beta"] = - y[,"S"]
    r[, "I", "S", "beta"] = y[, "I"]
    r[, "I", "I", "beta"] = y[, "S"]
    r[, "I", "I", "gamma"] = - 1
    return(r)
}

