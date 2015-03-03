ProfileSSE.covariance.delay <- function(g)
{
    H <- t(g)%*%g
    Covar <- NeweyWest.Var( 0.5*(t(H)+H) ,g,5)
    return( Covar )
}
