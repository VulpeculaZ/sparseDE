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

mDTVSIRtrfn <- list()

mDTVSIRtrfn$fn <- function (t, y, p, more)
{
    r = y
    yi.d <- more$y.d[,1]
    pk <- p[(length(p) - more$nKappa + 1):length(p)]
    phot <- p["pho0"] + t * p["pho1"]
    b <- more$b
    r[, "S"] =  phot - tvtrans(t, pk) * yi.d * y[, "S"] + b ## * p["alpha"]
    r[, "I"] =  tvtrans(t, pk) * yi.d * y[, "S"] - p["gamma"] * y[, "I"]
    return(r)
}

mDTVSIRtrfn$dfdx <- function (t, y, p, more)
{
    r = array(0, c(length(t), ncol(y), ncol(y)))
    dimnames(r) = list(NULL, colnames(y), colnames(y))
    pk <- p[(length(p) - more$nKappa + 1):length(p)]
    phot <- p["pho0"] + t * p["pho1"]
    yi.d <- more$y.d[,1]
    r[, "S", "S"] = - tvtrans(t, pk) * yi.d
    r[, "I", "S"] =  tvtrans(t, pk) * yi.d
    r[, "I", "I"] = -p["gamma"]
    return(r)
}

mDTVSIRtrfn$dfdx.d <- function (t, y, p, more)
{
    ## yi.d <- more$y.d[,2]
    pk <- p[(length(p) - more$nKappa + 1):length(p)]
    phot <- p["pho0"] + t * p["pho1"]
    r = array(0, c(length(t), ncol(y), ncol(y)))
    dimnames(r) = list(NULL, colnames(y), colnames(y))
    r[, "S", "I"] = - tvtrans(t, pk) * y[,"S"]
    r[, "I", "I"] =  tvtrans(t, pk) * y[,"S"]
    return(r)
}

mDTVSIRtrfn$dfdp <- function (t, y, p, more)
{
    ## b <- more$b
    phot <- p["pho0"] + t * p["pho1"]
    yi.d <- more$y.d[,1]
    r = array(0, c(length(t), ncol(y), length(p)))
    dimnames(r) = list(NULL, colnames(y), names(p))
    ## r[ , "S", "alpha"] = b
    r[, "I", "gamma"] = - y[, "I"]
    r[,"I", "pho0"] = 1
    r[,"I", "pho1"] = t
    month <- t %% 1
    nKappa <- more$nKappa
    r[ , "S", "k1"][month  >= (nKappa-1)/nKappa & month < 1] = -(yi.d * y[, "S"] * month)[month  >= (nKappa-1)/nKappa & month < 1]
    r[ , "S", "k1"][month  >= 0 & month < 1/nKappa] = -(yi.d * y[, "S"] * (1 - month))[month  >= 0 & month < 1/nKappa]
    r[ , "I", "k1"][month  >= (nKappa-1)/nKappa & month < 1] = (yi.d * y[, "S"] * month)[month  >= (nKappa-1)/nKappa & month < 1]
    r[ , "I", "k1"][month  >= (i-1)/nKappa & month < i/nKappa] = (yi.d * y[, "S"] * (1 - month))[month  >= 0 & month < 1/nKappa]
    for(i in 2 : nKappa){
        r[ , "S", paste("k", i, sep ="")][month  >= (i-2)/nKappa & month < (i-1)/nKappa] = -(yi.d * y[, "S"] * month)[month  >= (i-2)/nKappa & month < (i-1) /nKappa]
        r[ , "S", paste("k", i, sep ="")][month  >= (i-1)/nKappa & month < i/nKappa] = -(yi.d * y[, "S"] * (1 - month))[month  >= (i-1)/nKappa & month < i/nKappa]
        r[ , "I", paste("k", i, sep ="")][month  >= (i-2)/nKappa & month < (i-1)/nKappa] = (yi.d * y[, "S"] * month)[month  >= (i-2)/nKappa & month < (i-1) /nKappa]
        r[ , "I", paste("k", i, sep ="")][month  >= (i-1)/nKappa & month < i/nKappa] = (yi.d * y[, "S"] * (1 - month))[month  >= (i-1)/nKappa & month < i/nKappa]
    }
    return(r)
}

mDTVSIRtrfn$d2fdx2 <- function (t, y, p, more)
{
    r = array(0, c(length(t), ncol(y), ncol(y), ncol(y)))
    phot <- p["pho0"] + t * p["pho1"]
    dimnames(r) = list(NULL, colnames(y), colnames(y), colnames(y))
    return(r)
}

mDTVSIRtrfn$d2fdxdp <- function (t, y, p, more)
{
    yi.d <- more$y.d[,1]
    r = array(0, c(length(t), ncol(y), ncol(y), length(p)))
    phot <- p["pho0"] + t * p["pho1"]
    dimnames(r) = list(NULL, colnames(y), colnames(y), names(p))
    r[, "I", "I", "gamma"] = -1
    month <- t %% 1
    nKappa <- more$nKappa

    r[ , "S", "S", "k1"][month  >= (nKappa-1)/nKappa & month < 1] = -(yi.d * month)[month  >= (nKappa-1)/nKappa & month < 1]
    r[ , "S", "S","k1"][month  >= 0 & month < 1/nKappa] = -(yi.d * (1 - month))[month  >= 0 & month < 1/nKappa]
    r[ , "I", "S", "k1"][month  >= (nKappa-1)/nKappa & month < 1] = (yi.d * month)[month  >= (nKappa-1)/nKappa & month < 1]
    r[ , "I", "S", "k1"][month  >= (i-1)/nKappa & month < i/nKappa] = (yi.d * (1 - month))[month  >= 0 & month < 1/nKappa]
    for(i in 2 : nKappa){
        r[ , "S", "S", paste("k", i, sep ="")][month  >= (i-2)/nKappa & month < (i-1)/nKappa] = -(yi.d * month)[month  >= (i-2)/nKappa & month < (i-1) /nKappa]
        r[ , "S", "S", paste("k", i, sep ="")][month  >= (i-1)/nKappa & month < i/nKappa] = -(yi.d * (1 - month))[month  >= (i-1)/nKappa & month < i/nKappa]
        r[ , "I", "S", paste("k", i, sep ="")][month  >= (i-2)/nKappa & month < (i-1)/nKappa] = (yi.d * month)[month  >= (i-2)/nKappa & month < (i-1) /nKappa]
        r[ , "I", "S", paste("k", i, sep ="")][month  >= (i-1)/nKappa & month < i/nKappa] = (yi.d * (1 - month))[month  >= (i-1)/nKappa & month < i/nKappa]
    }
    return(r)
}


mDTVSIRtrfn$d2fdx.ddp <- function (t, y, p, more)
{
    yi.d <- more$y.d[,1]
    r = array(0, c(length(t), ncol(y), ncol(y), length(p)))
    phot <- p["pho0"] + t * p["pho1"]
    dimnames(r) = list(NULL, colnames(y), colnames(y), names(p))
    month <- t %% 1
    nKappa <- more$nKappa

    r[ , "S", "I", "k1"][month  >= (nKappa-1)/nKappa & month < 1] = -(y[, "S"] * month)[month  >= (nKappa-1)/nKappa & month < 1]
    r[ , "S", "I", "k1"][month  >= 0 & month < 1/nKappa] = -(y[, "S"] * (1 - month))[month  >= 0 & month < 1/nKappa]
    r[ , "I", "I", "k1"][month  >= (nKappa-1)/nKappa & month < 1] = (y[, "S"] * month)[month  >= (nKappa-1)/nKappa & month < 1]
    r[ , "I", "I", "k1"][month  >= (i-1)/nKappa & month < i/nKappa] = (y[, "S"] * (1 - month))[month  >= 0 & month < 1/nKappa]
    for(i in 2 : nKappa){
        r[ , "S", "I", paste("k", i, sep ="")][month  >= (i-2)/nKappa & month < (i-1)/nKappa] = -(y[, "S"] * month)[month  >= (i-2)/nKappa & month < (i-1) /nKappa]
        r[ , "S", "I", paste("k", i, sep ="")][month  >= (i-1)/nKappa & month < i/nKappa] = -(y[, "S"] * (1 - month))[month  >= (i-1)/nKappa & month < i/nKappa]
        r[ , "I", "I", paste("k", i, sep ="")][month  >= (i-2)/nKappa & month < (i-1)/nKappa] = (y[, "S"] * month)[month  >= (i-2)/nKappa & month < (i-1) /nKappa]
        r[ , "I", "I", paste("k", i, sep ="")][month  >= (i-1)/nKappa & month < i/nKappa] = (y[, "S"] * (1 - month))[month  >= (i-1)/nKappa & month < i/nKappa]
    }
    return(r)
}

mDTVSIRtrfn$d2fdxdx.d <- function (t, y, p, more)
{
    yi.d <- more$y.d[,1]
    pk <- p[(length(p) - more$nKappa + 1):length(p)]
    r = array(0, c(length(t), ncol(y), ncol(y), ncol(y)))
    dimnames(r) = list(NULL, colnames(y), colnames(y), colnames(y))
    r[, "S", "S", "I"] = - tvtrans(t, pk)
    r[, "I", "S", "I"] =  tvtrans(t, pk)
    return(r)
}


mDTVSIRtrfn$d2fdx.d2 <- function (t, y, p, more)
{
    yi.d <- more$y.d[,1]
    r = array(0, c(length(t), ncol(y), ncol(y), ncol(y)))
    phot <- p["pho0"] + t * p["pho1"]
    dimnames(r) = list(NULL, colnames(y), colnames(y), colnames(y))
    return(r)
}
