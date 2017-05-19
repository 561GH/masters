# MH within Gibbs for eta0

#' @import mvtnorm
#' @import stats
#' @export
logposterior <- function(l, sig2,
                          ystar, yprime,
                          S = NULL, Syy = NULL,
                          given) {

  # S: covariance matrix for (ystar, yprime); either provide or
  #   calculate accordingly
  # Syy: covariance matrix for y; either provide or
  #   calculate accordingly
  # given: list(y = y, xstar = xstar, xprime = xprime, x = x)
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

  if (is.null(S)) {
    S11 <- covMatrix(X = xstar, sig2 = sig2, covar.fun = r.matern, l = l)
    S22 <- covMatrix(X = xprime, sig2 = sig2, covar.fun = r.matern2, l = l)
    S21 <- covMatrix(X = xprime, X2 = xstar,  # CAREFUL SEE PAPER FOR ARG. ORDER
                     sig2 = sig2, covar.fun = r.matern1, l = l)
    S <- rbind(cbind(S11, t(S21)),
               cbind(S21, S22))

    # S21 <- covMatrix(X = xprime, X2 = xstar,
    #                  sig2 = sig2, covar.fun = r.matern1, l = l)
    # all.equal(S12, t(S21))  # should be TRUE...
  }

  if (is.null(Syy)) {
    Syy <- covMatrix(X = x, sig2 = sig2, covar.fun = r.matern, l = l)
  }

  ## use log scale because numbers are so tiny
  log.d.ystaryprime <- dmvnorm(c(ystar, yprime), mean = m, sigma = S, log = TRUE)
  log.d.y <- dmvnorm(as.numeric(y), mean = my, sigma = Syy, log = TRUE)

  # cat("log.d.ystaryprime:", log.d.ystaryprime, "\n")
  # cat("log.d.y:", log.d.y, "\n")

  return( log( dchisq(sqrt(5) / l, df = 1) ) + log( dchisq(sig2, df = 5) ) +
            log.d.ystaryprime + log.d.y )

}

