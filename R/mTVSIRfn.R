tvtrans <- function(t,k){
    month <- t %% 1
    r <- rep(0, length(t))
    for(i in 1:length(k)){
        r[month  >= (i-1)/length(k) & month < i /length(k)] <- k[i]
    }
    return(r)
}

mTVSIRfn <- list()

mTVSIRfn$fn <- function (t, y, p, more)
{
    r = y
    pk <- p[(length(p) - more$nKappa + 1):length(p)]
    b <- more$b
    r[, "S"] =  - tvtrans(t, pk) * y[,"I"] * y[, "S"] + b
    r[, "I"] =  tvtrans(t, pk) * y[,"I"] * y[, "S"] - p["gamma"] * y[, "I"]
    return(r)
}

mTVSIRfn$dfdx <- function (t, y, p, more)
{
    r = array(0, c(length(t), ncol(y), ncol(y)))
    dimnames(r) = list(NULL, colnames(y), colnames(y))
    pk <- p[(length(p) - more$nKappa + 1):length(p)]
    r[, "S", "S"] = - tvtrans(t, pk) * y[,"I"]
    r[, "I", "S"] =  tvtrans(t, pk) * y[,"I"]
    r[, "I", "I"] = -p["gamma"]
    return(r)
}

mTVSIRfn$dfdp <- function (t, y, p, more)
{
    r = array(0, c(length(t), ncol(y), length(p)))
    dimnames(r) = list(NULL, colnames(y), names(p))
    nKappa <- more$nKappa
    r[, "I", "gamma"] = - y[, "I"]
    b <- more$b
    month <- t %% 1
    for(i in 1:nKappa){
        r[ , "S", paste("k", i, sep ="")][month  >= (i-1)/nKappa & month < i /nKappa] = - (y[,"I"] * y[, "S"])[month  >= (i-1)/nKappa & month < i /nKappa]
        r[ , "I", paste("k", i, sep ="")][month  >= (i-1)/nKappa & month < i /nKappa] = (y[,"I"] * y[, "S"])[month  >= (i-1)/nKappa & month < i /nKappa]
    }
    return(r)
}

mTVSIRfn$d2fdx2 <- function (t, y, p, more)
{
    r = array(0, c(length(t), ncol(y), ncol(y), ncol(y)))
    dimnames(r) = list(NULL, colnames(y), colnames(y), colnames(y))
    return(r)
}

mTVSIRfn$d2fdxdp <- function (t, y, p, more)
{
    r = array(0, c(length(t), ncol(y), ncol(y), length(p)))
    dimnames(r) = list(NULL, colnames(y), colnames(y), names(p))
    nKappa <- more$nKappa
    r[, "I", "I", "gamma"] = -1
    month <- t %% 1
    for(i in 1:nKappa){
        r[ , "S", "S", paste("k", i, sep ="")][month  >= (i-1)/nKappa & month < i /nKappa] = - y[,"I"][month  >= (i-1)/nKappa & month < i /nKappa]
        r[ , "I", "S", paste("k", i, sep ="")][month  >= (i-1)/nKappa & month < i /nKappa] = y[,"I"][month  >= (i-1)/nKappa & month < i /nKappa]
        r[ , "S", "I", paste("k", i, sep ="")][month  >= (i-1)/nKappa & month < i /nKappa] = - y[,"S"][month  >= (i-1)/nKappa & month < i /nKappa]
        r[ , "I", "I", paste("k", i, sep ="")][month  >= (i-1)/nKappa & month < i /nKappa] = y[,"S"][month  >= (i-1)/nKappa & month < i /nKappa]
    }
    return(r)
}
