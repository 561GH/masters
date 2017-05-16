# MH within Gibbs for eta0

library(mvtnorm)

posterior <- function(l, sig2,
                      ystar, yprime, y,
                      xstar, xprime, x) {

  # xstar, xprime, x: careful likely have to cbind() themn to make them columns

  m <- c(rep(mean(ystar), length(ystar)),
         rep(mean(yprime), length(yprime)))
  my <- rep(mean(y), length(y))

  S <- covMatrix(X = rbind(xstar, xprime), sig2 = sig2,
                 covar.fun = r.matern, l = l)
  Syy <- covMatrix(X = x, sig2 = sig2,
                   covar.fun = r.matern, l = l)

  ## use log scale because numbers are so tiny
  log.d.ystaryprime <- dmvnorm(c(ystar, yprime), mean = m, sigma = S, log = TRUE)
  log.d.y <- dmvnorm(y, mean = my, sigma = Syy, log = TRUE)

  return( log( dchisq(sqrt(5) / l, df = 1) ) + log( dchisq(sig2, df = 5) ) +
            log.d.ystaryprime + log.d.y )

}
