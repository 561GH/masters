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

  nugget <- 1e-6
  R <- covMatrix(X1 = x, X2 = x, sig2 = sig2, l = l,
                 covar.fun = "matern") +
    diag(nugget, nrow(x))
  Linv <- solve( t(chol(R)) )  # recall need transpose to match std chol def
  Rinv <- t(Linv) %*% Linv

  S11 <- covMatrix(X1 = xstar, X2 = xstar,
                   sig2 = sig2, l = l, covar.fun = "matern")
  S22 <- covMatrix(X1 = xprime, X2 = xprime,
                   sig2 = sig2, l = l, covar.fun = "matern2", d1 = 1, d2 = 1)
  S21 <- covMatrix(X1 = xprime, X2 = xstar,  # CAREFUL SEE PAPER FOR ARG. ORDER
                   sig2 = sig2, l = l, covar.fun = "matern1", d1 = 1)

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
  # (I think the paper is wrong... following Rasmussen (2.19) instead)
  tau2.xstarprime <- R.xstarxprime -
    r.xstarprime.x %*% Rinv %*% t(r.xstarprime.x)

  # MAKE SURE tau2.xstarprime IS SYMMETRIC
  tau2.xstarprime <- ( tau2.xstarprime + t(tau2.xstarprime) ) / 2

  ## use log scale because numbers are so tiny
  # log.d.ystaryprime <- dmvnorm(c(ystar, yprime),
  #                              mean = mu.xstarprime, sigma = tau2.xstarprime,
  #                              log = TRUE)
  log.d.ystaryprime <- dmvn(c(ystar, yprime),
                               m = mu.xstarprime, S = tau2.xstarprime,
                               log = TRUE)
  log.d.y <- dmvnorm(as.numeric(y),
                     mean = rep(0, length(y)), sigma = R, log = TRUE)

  # cat("log.d.ystaryprime:", log.d.ystaryprime, "\n")
  # cat("log.d.y:", log.d.y, "\n")

  return( log( dchisq(sqrt(5) / l, df = 1) ) + log( dchisq(sig2, df = 5) ) +
            log.d.ystaryprime + log.d.y )

}

