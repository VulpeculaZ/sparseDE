DSIRfn.sparse <- list()

DSIRfn.sparse$fn <- function (t, y, p, more)
{
    r = y
    yi.d <- more$y.d[,1]
    r[, "S"] =  - yi.d * y[, "S"] + (1 - y[,"S"]) * (0.25 * (sin(t) + 1) + 0.1)
    r[, "I"] =  yi.d * y[, "S"] - (p["gamma"] + (0.25 * (sin(t) + 1) + 0.1)) * y[, "I"]
    return(r)
}

DSIRfn.sparse$dfdx <- function (t, y, p, more)
{
    r = array(0, c(length(t), ncol(y), ncol(y)))
    dimnames(r) = list(NULL, colnames(y), colnames(y))
    yi.d <- more$y.d[,1]
    r[, "S", "S"] = -  yi.d - (0.25 * (sin(t) + 1) + 0.1)
    r[, "I", "S"] = yi.d
    r[, "I", "I"] = -p["gamma"] - (0.25 * (sin(t) + 1) + 0.1)
    return(r)
}

DSIRfn.sparse$dfdx.d <- function (t, y, p, more)
{
    ## yi.d <- more$y.d[,2]
    r = array(0, c(length(t), ncol(y), ncol(y)))
    dimnames(r) = list(NULL, colnames(y), colnames(y))
    r[, "S", "I"] = - y[,"S"]
    r[, "I", "I"] = y[,"S"]
    return(r)
}



DSIRfn.sparse$dfdp <- function (t, y, p, more)
{
    yi.d <- more$y.d[,1]
    r = array(0, c(length(t), ncol(y), length(p)))
    dimnames(r) = list(NULL, colnames(y), names(p))
    r[, "I", "gamma"] = - y[, "I"]
    return(r)
}

DSIRfn.sparse$d2fdx2 <- function (t, y, p, more)
{
    r = array(0, c(length(t), ncol(y), ncol(y), ncol(y)))
    dimnames(r) = list(NULL, colnames(y), colnames(y), colnames(y))
    return(r)
}

DSIRfn.sparse$d2fdxdp <- function (t, y, p, more)
{
    yi.d <- more$y.d[,1]
    r = array(0, c(length(t), ncol(y), ncol(y), length(p)))
    dimnames(r) = list(NULL, colnames(y), colnames(y), names(p))
    r[, "I", "I", "gamma"] = - 1
    return(r)
}


DSIRfn.sparse$d2fdx.ddp <- function (t, y, p, more)
{
    yi.d <- more$y.d[,1]
    r = array(0, c(length(t), ncol(y), ncol(y), length(p)))
    dimnames(r) = list(NULL, colnames(y), colnames(y), names(p))
    return(r)
}

DSIRfn.sparse$d2fdxdx.d <- function (t, y, p, more)
{
    yi.d <- more$y.d[,1]
    r = array(0, c(length(t), ncol(y), ncol(y), ncol(y)))
    dimnames(r) = list(NULL, colnames(y), colnames(y), colnames(y))
    r[, "S", "S", "I"] = - 1
    r[, "I", "S", "I"] = 1
    return(r)
}


DSIRfn.sparse$d2fdx.d2 <- function (t, y, p, more)
{
    yi.d <- more$y.d[,1]
    r = array(0, c(length(t), ncol(y), ncol(y), ncol(y)))
    dimnames(r) = list(NULL, colnames(y), colnames(y), colnames(y))
    return(r)
}

