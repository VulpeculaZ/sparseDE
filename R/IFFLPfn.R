tvI <- function(t,k){
    rt <- c(t[1], t[length(t)])
    stepLen <- (rt[2] -rt[1]) / length(k)
    r <- rep(NA, length(t))
    for(i in 1 : length(k)){
        r[t  >= rt[1] + (i-1) * stepLen & t < rt[1] + i * stepLen] <- k[i]
    }
    r[length(t)] <- k[length(k)]
    return(r)
}

Hill <- function(L, K){
    L / (L + K)
}

dHdL <- function(L, K){
    K / (L + K)^2
}

dHdK <- function(L, K){
    (-1) * L / (L + K)^2
}

d2HdL2 <- function(L, K){
    (-2) * K / (L + K)^3
}

d2HdLdK <- function(L, K){
    (L - K) / (L + K)^3
}

IFFLPfn <- list()

IFFLPfn$fn <- function(t, y, p, more){
    r <- y
    pk <- p[(length(p) - more$nKappa + 1):length(p)]
    r[ , "A"] <- tvI(t, pk)  * Hill(1 - y[, "A"], p["KIA"]) - p["kFA"] * Hill(y[,"A"], p["KFA"])
    r[ , "B"] <- y[,"A"] * p["kAB"] * Hill(1 - y[, "B"], p["KAB"]) - p["kFB"] * Hill(y[,"B"], p["KFB"])
    r[ , "C"] <- y[,"A"] * p["kAC"] * Hill(1 - y[, "C"], p["KAC"]) - y["B"] * p["kBC"] * Hill(y[,"C"], p["KBC"])
    return(r)
}

IFFLPfn$dfdx <- function(t, y, p, more){
    r <- array(0, c(length(t), ncol(y), ncol(y)))
    dimnames(r) <- list(NULL, colnames(y), colnames(y))
    pk <- p[(length(p) - more$nKappa + 1):length(p)]
    r[ , "A", "A"] <- tvI(t, pk) * (-1) * dHdL(1 - y[, "A"], p["KIA"]) - p["kFA"] * dHdL(y[,"A"], p["KFA"])
    r[ , "B", "A"] <- p["kAB"] * Hill(1 - y[, "B"], p["KAB"])
    r[ , "B", "B"] <- y[,"A"] * p["kAB"] * (-1) * dHdL(1 - y[, "B"], p["KAB"]) - p["kFB"] * dHdL(y[,"B"], p["KFB"])
    r[ , "C", "A"] <- p["kAC"] * Hill(1 - y[, "C"], p["KAC"])
    r[ , "C", "B"] <- (-1) * p["kBC"] * Hill(y[,"C"], p["KBC"])
    r[ , "C", "C"] <- (-1) *  y[,"A"] * p["kAC"] * dHdL(1 - y[, "C"], p["KAC"]) - y["B"] * p["kBC"] * dHdL(y[,"C"], p["KBC"])
    return(r)
}

IFFLPfn$dfdp <- function(t, y, p, more){
    r = array(0, c(length(t), ncol(y), length(p)))
    dimnames(r) = list(NULL, colnames(y), names(p))
    nKappa <- more$nKappa
    pk <- p[(length(p) - nKappa + 1):length(p)]
    stepLen <- (t[length(t)] - t[1]) / nKappa

    r[ , "A", "kFA"] <- (-1) * Hill(y[,"A"], p["KFA"])
    r[ , "A", "KIA"] <-  p["kIA"] * tvI(t, pk) * dHdK(1 - y[, "A"], p["KIA"])
    r[ , "A", "KFA"] <- (-1) * p["kFA"] * dHdK(y[,"A"], p["KFA"])

    r[ , "B", "kAB"] <-  y[,"A"] * Hill(1 - y[, "B"], p["KAB"])
    r[ , "B", "kFB"] <- (-1) * Hill(y[,"B"], p["KFB"])
    r[ , "B", "KAB"] <-  p["kAB"] * y["A"] * dHdK(1 - y[, "B"], p["KAB"])
    r[ , "B", "KFB"] <- (-1) * p["kFB"] * dHdK(y[,"B"], p["KFB"])

    r[ , "C", "kAC"] <-  y[,"A"] * Hill(1 - y[, "C"], p["KAC"])
    r[ , "C", "kBC"] <- (-1) * y[ , "B"] * Hill(y[,"C"], p["KBC"])
    r[ , "C", "KAC"] <-  p["kAC"] * y["A"] * dHdK(1 - y[, "C"], p["KAC"])
    r[ , "C", "KBC"] <- (-1) * y["B"] * p["kBC"] * dHdK(y[,"C"], p["KBC"])

    dAdKappa <-  Hill(1 - y[, "A"], p["KIA"])
    for(i in 1:nKappa){
        r[ , "A", paste("k", i, sep ="")][t >= t[1] + (i-1) * stepLen & t < t[1] + i * stepLen] <-  dAdKappa[t >= t[1] + (i-1) * stepLen & t < t[1] + i * stepLen]
    }
    r[ length(t),"A",  paste("k", i, sep ="")] <- dAdKappa[length(t)]
    return(r)
}

IFFLPfn$d2fdx2 <- function(t, y, p, more){
    r = array(0, c(length(t), ncol(y), ncol(y), ncol(y)))
    dimnames(r) = list(NULL, colnames(y), colnames(y), colnames(y))
    pk <- p[(length(p) - more$nKappa + 1):length(p)]
    r[ , "A", "A"] <- tvI(t, pk)  * d2HdL2(1 - y[, "A"], p["KIA"]) - p["kFA"] * d2HdL2(y[,"A"], p["KFA"])
    r[ , "B", "A", "B"] <- r[ , "B", "B", "A"] <- p["kAB"] * (-1) * dHdL(1 - y[, "B"], p["KAB"])
    r[ , "B", "B", "B"] <- y[,"A"] * p["kAB"] * d2HdL2(1 - y[, "B"], p["KAB"]) - p["kFB"] * d2HdL2(y[,"B"], p["KFB"])
    r[ , "C", "A", "C"] <- r[ , "C", "C", "A"] <- (-1) * p["kAC"] *  dHdL(1 - y[, "C"], p["KAC"])
    r[ , "C", "B", "C"] <- r[ , "C", "C", "B"] <- (-1) * p["kBC"] * dHdL(y[,"C"], p["KBC"])
    r[ , "C", "C"] <- y[,"A"] * p["kAC"] * d2HdL2(1 - y[, "C"], p["KAC"]) - y["B"] * p["kBC"] * d2HdL2(y[,"C"], p["KBC"])
    return(r)
}

IFFLPfn$d2fdxdp <- function(t, y, p, more){
    r = array(0, c(length(t), ncol(y), ncol(y), length(p)))
    dimnames(r) = list(NULL, colnames(y), colnames(y), names(p))
    nKappa <- more$nKappa
    pk <- p[(length(p) - more$nKappa + 1):length(p)]
    stepLen <- (t[length(t)] - t[1]) / nKappa

    r[ , "A", "A", "kFA"] <- (-1) * dHdL(y[,"A"], p["KFA"])
    r[ , "A", "A", "KIA"] <- (-1) * tvI(t, pk) * d2HdLdK(1 - y[, "A"], p["KIA"])
    r[ , "A", "A", "KFA"] <- (-1) * p["kFA"] * d2HdLdK(y[,"A"], p["KFA"])

    r[ , "B", "A", "kAB"] <- Hill(1 - y[, "B"], p["KAB"])
    r[ , "B", "b", "kAB"] <- y["A"] * dHdL(1 - y[, "B"], p["KAB"]) * (-1)
    r[ , "B", "B", "kFB"] <- (-1) * dHdL(y[,"B"], p["KFB"])
    r[ , "B", "A", "KAB"] <- p["kAB"] * dHdK(1 - y[, "B"], p["KAB"])
    r[ , "B", "B", "KAB"] <- p["kAB"] * y["A"] * d2HdLdK(1 - y[, "B"], p["KAB"]) * (-1)
    r[ , "B", "B", "KFB"] <- (-1) * p["kFB"] * d2HdLdK(y[,"B"], p["KFB"])

    r[ , "C", "A", "kAC"] <- Hill(1 - y[, "C"], p["KAC"])
    r[ , "C", "C", "kAC"] <- y[,"A"] * dHdL(1 - y[, "C"], p["KAC"]) * (-1)
    r[ , "C", "B", "kBC"] <- (-1) * Hill(y[,"C"], p["KBC"])
    r[ , "C", "C", "kBC"] <- (-1) * y[ , "B"] * dHdL(y[,"C"], p["KBC"])

    r[ , "C", "A", "KAC"] <- p["kAC"] * dHdK(1 - y[, "C"], p["KAC"])
    r[ , "C", "C", "KAC"] <-  p["kAC"] * y["A"] * d2HdLdK(1 - y[, "C"], p["KAC"]) * (-1)
    r[ , "C", "B", "KBC"] <- (-1) * p["kBC"] * dHdK(y[,"C"], p["KBC"])
    r[ , "C", "C", "KBC"] <- (-1) * y["B"] * p["kBC"] * d2HdLdK(y[,"C"], p["KBC"])

    d2AdKappadA <- (-1) * dHdL(1 - y[, "A"], p["KIA"])
    for(i in 1:nKappa){
        r[ , "A", "A", paste("k", i, sep ="")][t >= t[1] + (i-1) * stepLen & t < t[1] + i * stepLen] <-  dAdKappa[t >= t[1] + (i-1) * stepLen & t < t[1] + i * stepLen]
    }
    r[ length(t),"A", "A", paste("k", i, sep ="")] <- dAdKappa[length(t)]
    return(r)
}



