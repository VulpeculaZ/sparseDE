##' .. content for \description{} (no empty lines) ..
##' A function to plot time-series diagnostics of the residuals from a generalized profiling DDE model.
##' .. content for \details{} ..
##' @title Time Series Diagnostics for the Residuals
##' @param y Matrix of observed data values.
##' @param times Vector observation times for the data.
##' @param fitted The functional data object for the estimated state process.
##' @param use.TSA \code{TURE} or \code{FALSE}. Whether to use TSA package to perfor the diagnostics.
##' @param ... Additional arguments for \code{tsdiag} function.
##' @return None. Diagnostics are plotted.
##' @seealso \code{tsdiag}
##' @author Ziqian Zhou
DDEdiag <- function(y, times, fitted, use.TSA = FALSE, ...){
    res <- y - eval.fd(times.d, fitted)
    if ((ncol(y)) > 1)
        cat("Multiple plots:  Click in the plot to advance to the next plot\n")

    for(i in 1:ncol(y)){
        if(use.TSA == FALSE){
            res.arima <- arima(res[,i], include.mean = FALSE)
            tsdiag(res.arima, ...)
        }
        else{
            res.arima <- arima(res[,i], include.mean = FALSE)
            tsdiag.Arima(res.arima, ...)
        }
    }
}
