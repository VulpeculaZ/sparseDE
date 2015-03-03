blowfliesfn <- list()

blowfliesfn$fn <- function(t, y, p, more){
    r = y
    y.d <- more$y.d[,1]
    r[,"y"] <- p["c"] * y.d * exp(-y.d / p["N0"]) - p["a"] * y
    return(r)
}

blowfliesfn$dfdx <- function(t, y, p, more){
    r <- array(0, c(length(t), ncol(y), ncol(y)))
    dimnames(r) = list(NULL, colnames(y), colnames(y))
    r[,"y", "y"] <- -p["a"]
    return(r)
}

blowfliesfn$dfdx.d <- function(t, y, p, more){
    r <- array(0, c(length(t), ncol(y), ncol(y)))
    dimnames(r) = list(NULL, colnames(y), colnames(y))
    y.d <- more$y.d[,1]
    r[,"y", "y"] <- p["c"] * (1 - y.d / p["N0"]) * exp(-y.d / p["N0"])
    return(r)
}

blowfliesfn$dfdp <- function(t, y, p, more){
    r <- array(0, c(length(t), ncol(y), length(p)))
    dimnames(r) = list(NULL, colnames(y), names(p))
    y.d <- more$y.d[,1]
    r[,"y", "c"] <- y.d * exp(-y.d / p["N0"])
    r[, "y", "N0"] <- p["c"] * (y.d / p["N0"])^2 * exp(-y.d / p["N0"])
    r[, "y", "a"] <- -y
    return(r)
}

blowfliesfn$d2fdx2 <- function(t, y, p, more){
    r <- array(0, c(length(t), ncol(y), ncol(y), ncol(y)))
    dimnames(r) = list(NULL, colnames(y), colnames(y), colnames(y))
    return(r)
}

blowfliesfn$d2fdxdp <- function(t, y, p, more){
    r <- array(0, c(length(t), ncol(y), ncol(y), length(p)))
    dimnames(r) = list(NULL, colnames(y), colnames(y), names(p))
    y.d <- more$y.d[,1]
    r[,"y", "y", "a"] <- -1
    return(r)
}

blowfliesfn$d2fdx.ddp <- function(t, y, p, more){
    r <- array(0, c(length(t), ncol(y), ncol(y), length(p)))
    dimnames(r) = list(NULL, colnames(y), colnames(y), names(p))
    y.d <- more$y.d[,1]
    r[,"y", "y", "c"] <- exp(-y.d / p["N0"]) * (1 - y.d / p["N0"])
    r[,"y", "y", "N0"] <- p["c"] * exp(-y.d / p["N0"]) * (2 * y.d / p["N0"]^2 - y.d^2 / p["N0"]^3)
    return(r)
}

blowfliesfn$d2fdx.d2 <- function(t, y, p, more){
    r <- array(0, c(length(t), ncol(y), ncol(y), ncol(y)))
    dimnames(r) = list(NULL, colnames(y), colnames(y), colnames(y))
    y.d <- more$y.d[,1]
    r[,"y","y","y"] <- p["c"] * exp(-y.d / p["N0"]) * (y.d / p["N0"]^2 - 2 / p["N0"])
    return(r)
}

blowfliesfn$d2fdxdx.d <- function(t, y, p, more){
    r <- array(0, c(length(t), ncol(y), ncol(y), ncol(y)))
    dimnames(r) = list(NULL, colnames(y), colnames(y), colnames(y))
    return(r)
}

