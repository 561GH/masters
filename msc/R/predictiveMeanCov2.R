# Conditional mean and variance calculations for dim(x) = 2
# i.e. mean and variance of (y*, yprime) | y, l, sig2
#' @export
predictiveMeanCov2 <- function(given, l, sig2, nugget1 = 0, nugget2 = 0) {
  # given: x, xstar, xprime, y; y must be drawn from 0 mean GP
  # l: length-scale GP parameter
  # sig2: constant variance GP parameter

  x <- given$x
  xstar <- given$xstar
  xprime <- given$xprime
  y <- given$y

  R <- covMatrix(X1 = x, X2 = x, sig2 = sig2, l = l,
                 covar.fun = "matern") +
    diag(nugget1, nrow(x))
  Linv <- solve( t(chol(R)) )  # recall need transpose to match std chol def
  Rinv <- t(Linv) %*% Linv

  n1 <- nrow(xprime) / 2  # number of derivs. in d1
  xprime1 <- xprime[1:n1,]
  xprime2 <- xprime[-(1:n1),]

  S11 <- covMatrix(X1 = xstar, X2 = xstar,
                   sig2 = sig2, l = l, covar.fun = "matern")
  S22 <- covMatrix(X1 = xprime1, X2 = xprime1,
                   sig2 = sig2, l = l, covar.fun = "matern2", d1 = 1, d2 = 1)
  S33 <- covMatrix(X1 = xprime2, X2 = xprime2,
                   sig2 = sig2, l = l, covar.fun = "matern2", d1 = 2, d2 = 2)

  S21 <- covMatrix(X1 = xprime1, X2 = xstar, # CAREFUL SEE PAPER FOR ARG. ORDER
                   sig2 = sig2, l = l, covar.fun = "matern1", d1 = 1)
  S31 <- covMatrix(X1 = xprime2, X2 = xstar, # CAREFUL SEE PAPER FOR ARG. ORDER
                   sig2 = sig2, l = l, covar.fun = "matern1", d1 = 2)

  S32 <- covMatrix(X1 = xprime2, X2 = xprime1,
                   sig2 = sig2, l = l, covar.fun = "matern2", d1 = 2, d2 = 1)

  R.xstarxprime <- rbind(cbind(S11, t(S21), t(S31)),
                         cbind(S21,   S22,  t(S32)),
                         cbind(S31,   S32,    S33)) +
    diag(rep(nugget2, nrow(xstar) + nrow(xprime)))

  # CAREFUL SEE PAPER FOR ARG. ORDER
  r.xstarprime.x <- rbind(covMatrix(X1 = xstar, X2 = x,
                                    sig2 = sig2, l = l,
                                    covar.fun = "matern"),
                          covMatrix(X1 = xprime1, X2 = x,
                                    sig2 = sig2, l = l,
                                    covar.fun = "matern1", d1 = 1),
                          covMatrix(X1 = xprime2, X2 = x,
                                    sig2 = sig2, l = l,
                                    covar.fun = "matern1", d1 = 2))
  mu.xstarprime <- r.xstarprime.x %*% Rinv %*% y

  # Following Rasmussen (2.19)
  tau2.xstarprime <- R.xstarxprime -
    r.xstarprime.x %*% Rinv %*% t(r.xstarprime.x) +
    diag(rep(0, nrow(xstar) + nrow(xprime)))

  # MAKE SURE tau2.xstarprime IS SYMMETRIC
  tau2.xstarprime <- ( tau2.xstarprime + t(tau2.xstarprime) ) / 2


  return(list(mu.xstarprime = mu.xstarprime,
              cov.xstarprime = tau2.xstarprime,
              R.xx = R))
}
