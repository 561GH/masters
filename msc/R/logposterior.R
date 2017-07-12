# MH within Gibbs for eta0

#' @import mvtnorm
#' @import stats
#' @export
logposterior <- function(l,
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
  m <- rep(0, length(ystar) + length(yprime))
  my <- rep(0, length(y)) #rep(mean(y), length(y))

  S11 <- covMatrix(X1 = xstar, X2 = xstar,
                   sig2 = sig2, l = l, covar.fun = "matern")
  S22 <- covMatrix(X1 = xprime, X2 = xprime,
                   sig2 = sig2, l = l, covar.fun = "matern2", d1 = 1, d2 = 1)
  S21 <- covMatrix(X1 = xprime, X2 = xstar,  # CAREFUL SEE PAPER FOR ARG. ORDER
                   sig2 = sig2, l = l, covar.fun = "matern1", d1 = 1)
  S <- rbind(cbind(S11, t(S21)),
             cbind(S21, S22))

  Syy <- covMatrix(X1 = x, X2 = x,
                   sig2 = sig2, l = l, covar.fun = "matern")

  ## use log scale because numbers are so tiny
  log.d.ystaryprime <- dmvnorm(c(ystar, yprime), mean = m, sigma = S, log = TRUE)
  log.d.y <- dmvnorm(as.numeric(y), mean = my, sigma = Syy, log = TRUE)

  # cat("log.d.ystaryprime:", log.d.ystaryprime, "\n")
  # cat("log.d.y:", log.d.y, "\n")

  return( log( dchisq(sqrt(5) / l, df = 1) ) + log( dchisq(sig2, df = 5) ) +
            log.d.ystaryprime + log.d.y )

}

