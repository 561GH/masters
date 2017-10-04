# Conditional mean and variance calculations
# i.e. mean and variance of (y*, yprime) | y, l, sig2
#' @export
predictiveMeanCov <- function(given, l, sig2) {
  # given: x, xstar, xprime, y; y must be drawn from 0 mean GP
  # l: length-scale GP parameter
  # sig2: constant variance GP parameter

  x <- given$x
  xstar <- given$xstar
  xprime <- given$xprime
  y <- given$y

  nugget <- 0#1e-6
  R <- covMatrix(X1 = x, X2 = x, sig2 = sig2, l = l,
                 covar.fun = "matern") +
    diag(nugget, nrow(x))
  Linv <- solve( t(chol(R)) )  # recall need transpose to match std chol def
  Rinv <- t(Linv) %*% Linv

  S11 <- covMatrix(X1 = xstar, X2 = xstar, sig2 = sig2, l = l,
                   covar.fun = "matern")
  S22 <- covMatrix(X1 = xprime, X2 = xprime, sig2 = sig2, l = l,
                   covar.fun = "matern2", d1 = 1, d2 = 1)
  S21 <- covMatrix(X1 = xprime, X2 = xstar,  # CAREFUL SEE PAPER FOR ARG. ORDER
                   sig2 = sig2, l = l,
                   covar.fun = "matern1", d1 = 1)
  R.xstarxprime <- rbind(cbind(S11, t(S21)),
                         cbind(S21, S22)) +
    diag(nugget, nrow(xstar) + nrow(xprime))

  # CAREFUL SEE PAPER FOR ARG. ORDER
  r.xstarprime.x <- rbind(covMatrix(X1 = xstar, X2 = x,
                                    sig2 = sig2, l = l,
                                    covar.fun = "matern"),
                          covMatrix(X1 = xprime, X2 = x,
                                    sig2 = sig2, l = l,
                                    covar.fun = "matern1", d1 = 1))
  mu.xstarprime <- r.xstarprime.x %*% Rinv %*% y

  # because covMatrix has sig2 in it, have an extra one multiplying in???
  tau2.xstarprime <- R.xstarxprime -
    r.xstarprime.x %*% Rinv %*% t(r.xstarprime.x) #/ sig2current #???

  # MAKE SURE tau2.xstarprime IS SYMMETRIC
  tau2.xstarprime <- ( tau2.xstarprime + t(tau2.xstarprime) ) / 2

  return(list(mu.xstarprime = mu.xstarprime,
              cov.xstarprime = tau2.xstarprime,
              R.xx = R))
}
