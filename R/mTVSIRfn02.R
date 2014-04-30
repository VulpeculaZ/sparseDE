## Continuous and ensured tvtrans(0,k) == tvtrans(1,k)
tvtrans <- function(t,k){
    month <- t %% 1
    r <- rep(0, length(t))
    ka <- c(k, k[1])
    for(i in 1:length(k)){
        mk <- month[month  >= (i-1)/length(k) & month < i /length(k)]
        r[month  >= (i-1)/length(k) & month < i /length(k)] <- k[i] + (ka[i+1] - k[i]) * (mk - (i-1)/length(k)) * length(k)
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
    nKappa <- more$nKappa
    dimnames(r) = list(NULL, colnames(y), colnames(y))
    pk <- p[(length(p) - nKappa + 1):length(p)]
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
    pk <- p[(length(p) - nKappa + 1):length(p)]
    month <- t %% 1
    for(i in 1:nKappa){
        mk <- month[month  >= (i-1)/length(k) & month < i /length(k)]
        r[ , "S", paste("k", i, sep ="")][month  >= (i-1)/nKappa & month < i /nKappa] = - (y[,"I"] * y[, "S"])[month  >= (i-1)/nKappa & month < i /nKappa] * (1- (mk - (i-1)/nKappa) * nKappa)
        r[ , "S", paste("k", (i+1) %% nKappa, sep ="")][month  >= (i-1)/nKappa & month < i /nKappa] <- - (y[,"I"] * y[, "S"])[month  >= (i-1)/nKappa & month < i /nKappa] * (mk - (i-1)/nKappa) * nKappa
        r[ , "I", paste("k", i, sep ="")][month  >= (i-1)/nKappa & month < i /nKappa] = (y[,"I"] * y[, "S"])[month  >= (i-1)/nKappa & month < i /nKappa] * (1- (mk - (i-1)/nKappa) * nKappa)
        r[ , "I", paste("k", (i+1) %% nKappa, sep ="")][month  >= (i-1)/nKappa & month < i /nKappa] <- (y[,"I"] * y[, "S"])[month  >= (i-1)/nKappa & month < i /nKappa] * (mk - (i-1)/nKappa) * nKappa

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
        mk <- month[month  >= (i-1)/length(k) & month < i /length(k)]
        r[ , "S"," S", paste("k", i, sep ="")][month  >= (i-1)/nKappa & month < i /nKappa] = - (y[,"I"]])[month  >= (i-1)/nKappa & month < i /nKappa] * (1- (mk - (i-1)/nKappa) * nKappa)
        r[ , "S", "S", paste("k", (i+1) %% nKappa, sep ="")][month  >= (i-1)/nKappa & month < i /nKappa] <- - (y[,"I"])[month  >= (i-1)/nKappa & month < i /nKappa] * (mk - (i-1)/nKappa) * nKappa
        r[ , "I", "S", paste("k", i, sep ="")][month  >= (i-1)/nKappa & month < i /nKappa] = (y[,"I"])[month  >= (i-1)/nKappa & month < i /nKappa] * (1- (mk - (i-1)/nKappa) * nKappa)
        r[ , "I", "S",  paste("k", (i+1) %% nKappa, sep ="")][month  >= (i-1)/nKappa & month < i /nKappa] <- (y[,"I"])[month  >= (i-1)/nKappa & month < i /nKappa] * (mk - (i-1)/nKappa) * nKappa
        r[ , "S", "I", paste("k", i, sep ="")][month  >= (i-1)/nKappa & month < i /nKappa] = - (y[, "S"])[month  >= (i-1)/nKappa & month < i /nKappa] * (1- (mk - (i-1)/nKappa) * nKappa)
        r[ , "S","I", paste("k", (i+1) %% nKappa, sep ="")][month  >= (i-1)/nKappa & month < i /nKappa] <- - (y[, "S"])[month  >= (i-1)/nKappa & month < i /nKappa] * (mk - (i-1)/nKappa) * nKappa
        r[ , "I","I", paste("k", i, sep ="")][month  >= (i-1)/nKappa & month < i /nKappa] = (y[, "S"])[month  >= (i-1)/nKappa & month < i /nKappa] * (1- (mk - (i-1)/nKappa) * nKappa)
        r[ , "I","I", paste("k", (i+1) %% nKappa, sep ="")][month  >= (i-1)/nKappa & month < i /nKappa] <- (y[, "S"])[month  >= (i-1)/nKappa & month < i /nKappa] * (mk - (i-1)/nKappa) * nKappa
    }
    return(r)
}
