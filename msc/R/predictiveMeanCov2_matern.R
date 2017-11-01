# Conditional mean and variance calculations for dim(x) = 2
# i.e. mean and variance of (y*, yprime) | y, l, sig2
#' @export
predictiveMeanCov2 <- function(given, l, sig2) {
  # given: x, xstar, xprime, y; y must be drawn from 0 mean GP
  # l: length-scale GP parameter
  # sig2: constant variance GP parameter

  x <- given$x
  xstar <- given$xstar
  xprime <- given$xprime
  y <- given$y

  # separate xprime for each dimension where want to take derivatives
  n1 <- nrow(xprime) / 2  # number of derivs. in d1
  xprime1 <- xprime[1:n1,]
  xprime2 <- xprime[-(1:n1),]

  ## K is 4 x 4 blocks
  k11 <- covMatrix(X1 = xstar, X2 = xstar, sig2 = sig2, l = l,
                   covar.fun = "matern")
  k44 <- covMatrix(X1 = x, X2 = x, sig2 = sig2, l = l,
                   covar.fun = "matern")
  k41 <- covMatrix(X1 = x, X2 = xstar, sig2 = sig2, l = l,
                   covar.fun = "matern")

  k22 <- covMatrix(X1 = xprime1, X2 = xprime1, sig2 = sig2, l = l,
                   covar.fun = "matern2", d1 = 1, d2 = 1)
  k33 <- covMatrix(X1 = xprime2, X2 = xprime2, sig2 = sig2, l = l,
                   covar.fun = "matern2", d1 = 2, d2 = 2)
  k32 <- covMatrix(X1 = xprime2, X2 = xprime1, sig2 = sig2, l = l,
                   covar.fun = "matern2", d1 = 2, d2 = 1)


  k21 <- covMatrix(X1 = xprime1, X2 = xstar, sig2 = sig2, l = l,
                   covar.fun = "matern1", d1 = 1)
  k31 <- covMatrix(X1 = xprime2, X2 = xstar, sig2 = sig2, l = l,
                   covar.fun = "matern1", d1 = 2)
  k42 <- covMatrix(X1 = x, X2 = xprime1, sig2 = sig2, l = l,
                   covar.fun = "matern1", d2 = 1)
  k43 <- covMatrix(X1 = x, X2 = xprime2, sig2 = sig2, l = l,
                   covar.fun = "matern1", d2 = 2)

  # NOTE: can use transpose for upper triangular part of K
  # only when derivative input sets have > 1 point
  # i.e. none of the k entries are single numbers
  # if single numbers, need to take negative for some of them
  # (see morris_simple_example C matrix)
  #
  # K <- rbind( cbind(k11, t(k21), t(k31), t(k41)),
  #             cbind(k21,   k22 , t(k32), t(k42)),
  #             cbind(k31,   k32 ,   k33 , t(k43)),
  #             cbind(k41,   k42 ,   k43 ,   k44) )

  S11 <- rbind( cbind(k11, t(k21), t(k31)),
                cbind(k21,   k22 , t(k32)),
                cbind(k31,   k32 ,   k33 ) )
  S12 <- rbind(t(k41),
               t(k42),
               t(k43))
  S22 <- k44

  # S22 == R matrix in Golchi
  # Sig1 == tau2.xstarprime matrix in Golchi
  Linv <- solve( t(chol(S22)) )  # recall need transpose to match std chol def
  S22inv <- t(Linv) %*% Linv
  Sig1 <- S11 - S12 %*% S22inv %*% t(S12)
  # Sig1 <- S11 - S12 %*% solve(S22, t(S12))  # not using chol decomposition

  # recall we have a 0 mean GP
  m <- S12 %*% S22inv %*% y

  return(list(mu.xstarprime = m,
              cov.xstarprime = Sig1,
              R.xx = S22))
}
