# MH within Gibbs for eta0

library(mvtnorm)

log.posterior <- function(l, sig2,
                          ystar, yprime,
                          given) {

  # given = list(y = y, xstar = xstar, xprime = xprime, x = x)
  # xstar, xprime, x: careful likely have to cbind() themn to make them columns

  ## unpacking given (cbind to make sure they are matrices)
  y <- given$y
  xstar <- cbind(given$xstar)
  xprime <- cbind(given$xprime)
  x <- cbind(given$x)

  ## means and covariances
  m <- c(rep(mean(ystar), length(ystar)),
         rep(mean(yprime), length(yprime)))
  my <- rep(mean(y), length(y))

  S11 <- covMatrix(X = xstar, sig2 = sig2, covar.fun = r.matern, l = l)
  S22 <- covMatrix(X = xprime, sig2 = sig2, covar.fun = r.matern2, l = l)
  S12 <- covMatrix(X = xstar, X2 = xprime,
                   sig2 = sig2, covar.fun = r.matern1, l = l)
  S <- rbind(cbind(S11, S12),
             cbind(t(S12), S22))
  # S21 <- covMatrix(X = xprime, X2 = xstar,
  #                  sig2 = sig2, covar.fun = r.matern1, l = l)
  # all.equal(S12, t(S21))  # should be TRUE...

  Syy <- covMatrix(X = x, sig2 = sig2, covar.fun = r.matern, l = l)

  ## use log scale because numbers are so tiny
  log.d.ystaryprime <- dmvnorm(c(ystar, yprime), mean = m, sigma = S, log = TRUE)
  log.d.y <- dmvnorm(y, mean = my, sigma = Syy, log = TRUE)

  cat("log.d.ystaryprime:", log.d.ystaryprime)
  cat("log.d.y:", log.d.y)

  return( log( dchisq(sqrt(5) / l, df = 1) ) + log( dchisq(sig2, df = 5) ) +
            log.d.ystaryprime + log.d.y )

}

