
DTVSIRfn.make <- function(){
    fn <- function (t, y, p, more)
    {
        r = y
        yi.d <- more$y.d[,1]
        pk <- p[(length(p) - more$nKappa + 1):length(p)]
        r[, "S"] =  - tvtrans(t, pk) * yi.d * y[, "S"] + 4000 * (sin(t / pi) + 2)
        r[, "I"] =  tvtrans(t, pk) * yi.d * y[, "S"] - p["gamma"] * y[, "I"]
        return(r)
    }

    dfdx <- function (t, y, p, more)
    {
        r = array(0, c(length(t), ncol(y), ncol(y)))
        dimnames(r) = list(NULL, colnames(y), colnames(y))
        pk <- p[(length(p) - more$nKappa + 1):length(p)]
        yi.d <- more$y.d[,1]
        r[, "S", "S"] = - tvtrans(t, pk) * yi.d
        r[, "I", "S"] =  tvtrans(t, pk) * yi.d
        r[, "I", "I"] = -p["gamma"]
        return(r)
    }

    dfdx.d <- function (t, y, p, more)
    {
        pk <- p[(length(p) - more$nKappa + 1):length(p)]
        r = array(0, c(length(t), ncol(y), ncol(y)))
        dimnames(r) = list(NULL, colnames(y), colnames(y))
        r[, "S", "I"] = - tvtrans(t, pk) * y[,"S"]
        r[, "I", "I"] =  tvtrans(t, pk) * y[,"S"]
        return(r)
    }

    dfdp <- function (t, y, p, more)
    {
        yi.d <- more$y.d[,1]
        r = array(0, c(length(t), ncol(y), length(p)))
        dimnames(r) = list(NULL, colnames(y), names(p))
        r[, "I", "gamma"] = - y[, "I"]
        period <- t %% 1
        nKappa <- more$nKappa
        for(i in 1 : nKappa){
            r[ , "S", paste("k", i, sep ="")][period  >= (i-1)/nKappa & period < i /nKappa] =
                - (yi.d * y[, "S"])[period  >= (i-1)/nKappa & period < i /nKappa]
            r[ , "I", paste("k", i, sep ="")][period  >= (i-1)/nKappa & period < i /nKappa] =
                (yi.d * y[, "S"])[period  >= (i-1)/nKappa & period < i /nKappa]
        }
        return(r)
    }

    d2fdx2 <- function (t, y, p, more)
    {
        r = array(0, c(length(t), ncol(y), ncol(y), ncol(y)))
        dimnames(r) = list(NULL, colnames(y), colnames(y), colnames(y))
        return(r)
    }

    d2fdxdp <- function (t, y, p, more)
    {
        yi.d <- more$y.d[,1]
        r = array(0, c(length(t), ncol(y), ncol(y), length(p)))
        dimnames(r) = list(NULL, colnames(y), colnames(y), names(p))
        r[, "I", "I", "gamma"] = -1
        period <- (t %% 5) / 5
        nKappa <- more$nKappa
        for(i in 1:nKappa){
            r[ , "S", "S", paste("k", i, sep ="")][period  >= (i-1)/nKappa & period < i /nKappa] = - yi.d[period  >= (i-1)/nKappa & period < i /nKappa]
            r[ , "I", "S", paste("k", i, sep ="")][period  >= (i-1)/nKappa & period < i /nKappa] = yi.d[period  >= (i-1)/nKappa & period < i /nKappa]
        }
        return(r)
    }

    d2fdx.ddp <- function (t, y, p, more)
    {
        yi.d <- more$y.d[,1]
        r = array(0, c(length(t), ncol(y), ncol(y), length(p)))
        dimnames(r) = list(NULL, colnames(y), colnames(y), names(p))
        period <- (t %% 5) / 5
        nKappa <- more$nKappa
        for(i in 1:nKappa){
            r[ , "S", "I", paste("k", i, sep ="")][period  >= (i-1)/nKappa & period < i /nKappa] = - y[,"S"][period  >= (i-1)/nKappa & period < i /nKappa]
            r[ , "I", "I", paste("k", i, sep ="")][period  >= (i-1)/nKappa & period < i /nKappa] = y[,"S"][period  >= (i-1)/nKappa & period < i /nKappa]
        }
        return(r)
    }


    d2fdxdx.d <- function (t, y, p, more)
    {
        yi.d <- more$y.d[,1]
        pk <- p[(length(p) - more$nKappa + 1):length(p)]
        r = array(0, c(length(t), ncol(y), ncol(y), ncol(y)))
        dimnames(r) = list(NULL, colnames(y), colnames(y), colnames(y))
        r[, "S", "S", "I"] = - tvtrans(t, pk)
        r[, "I", "S", "I"] =  tvtrans(t, pk)
        return(r)
    }


    d2fdx.d2 <- function (t, y, p, more)
    {
        yi.d <- more$y.d[,1]
        r = array(0, c(length(t), ncol(y), ncol(y), ncol(y)))
        dimnames(r) = list(NULL, colnames(y), colnames(y), colnames(y))
        return(r)
    }

    return(list(
        fn = fn, dfdx = dfdx,
        dfdp = dfdp, d2fdx2 = d2fdx2,
        d2fdxdp = d2fdxdp, d2fdx.ddp = d2fdx.ddp,
        dfdx.d = dfdx.d, d2fdx.ddx = d2fdxdx.d,
        d2fdxdx.d = d2fdxdx.d, d2fdx.d2 = d2fdx.d2
        ))
}
