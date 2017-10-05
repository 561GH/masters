# MH within Gibbs for eta0

#' @import stats
#' @export
logposterior2 <- function(l,
                         sig2,
                         ystar,
                         yprime,
                         given) {

  # given: list(y = y, xstar = xstar, xprime = xprime, x = x)
  # xstar, xprime, x: careful likely have to cbind() themn to make them columns

  ## unpacking given (cbind to make sure they are matrices)
  y <- given$y
  xstar <- cbind(given$xstar)
  xprime <- cbind(given$xprime)
  x <- cbind(given$x)

  ## means and covariances
  nugget <- 0
  kriging <- predictiveMeanCov2(given = given, l = l, sig2 = sig2)

  mu.xstarprime <- kriging$mu.xstarprime
  tau2.xstarprime <- kriging$cov.xstarprime +
    diag(rep(2, nrow(xstar) + nrow(xprime)))

  ## use log scale because numbers are so tiny
  log.d.y <- dmvn(as.numeric(y),
                  m = rep(0, length(y)), S = kriging$R.xx,
                  log = TRUE)

  log.d.ystaryprime <- dmvn(c(ystar, yprime),
                            m = mu.xstarprime, S = tau2.xstarprime,
                            log = TRUE)

  return( sum( log( dchisq(sqrt(5) / l, df = 1) ) ) +
            log( dchisq(sig2, df = 5) ) +
            log.d.ystaryprime + log.d.y )

}

