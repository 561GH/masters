## run setup

# monotone increasing function
ytrue <- function(x) {log(20*x + 1)}  # x > -1/20

# inputs where fcn evaluated
x <- cbind(c(0, 0.1, 0.2, 0.3, 0.4, 0.9, 1))
# true/known model values
y <- ytrue(x)

# prediction set inputs
xstar <- cbind(seq(0, 1, length.out = 50))

# derivative set inputs
xprime <- cbind(c(0.42, 0.47, 0.52, 0.57, 0.62, 0.67, 0.72, 0.77, 0.82, 0.87))

given <- list(x = cbind(x), xprime = cbind(xprime), xstar = cbind(xstar),
              y = y)

## start MH within gibbs

eta.init <- list(l = 10, sig2 = 1,
                 ystar = ytrue(xstar),
                 yprime = 20 / (20 * xprime))

set.seed(761)
N <- 10
v1 <- 1
v2 <- 1
chain.l <- rep(NA, N)
chain.sig2 <- rep(NA, N)
chain.ystar <- matrix(NA, nrow = length(xstar), ncol = N)
chain.yprime <- matrix(NA, nrow = length(xprime), ncol = N)

#out <- capture.output(
for (i in 1:N) {

  cat("Starting iteration:", i, "\n")

  if (i == 1) {
    lold <- eta.init$l
    sig2old <- eta.init$sig2
    ystarold <- eta.init$ystar
    yprimeold <- eta.init$yprime
  } else {
    lold <- chain.l[i-1]
    sig2old <- chain.l[i-1]
    ystarold <- chain.ystar[,i-1]
    yprimeold <- chain.yprime[,i-1]
  }

  ## UPDATE l ##################################################################
  lnew <- rnorm(1, mean = lold, sd = v1)

  if (lnew > 0) {

    logHR <- logposterior(l = lnew,
                           sig2 = sig2old, ystar = ystarold, yprime = yprimeold,
                           given = given) -
      logposterior(l = lold,
                    sig2 = sig2old, ystar = ystarold, yprime = yprimeold,
                    given = given)

    cat("lnew logHR:", logHR, "\n")

    if ( !is.nan(logHR) & !is.na(logHR) ) {  # when have -Inf - Inf get NaN

      if ( logHR > log(runif(1)) )  {
        chain.l[i] <- lnew
      } else {
        chain.l[i] <- lold
      }

    }

  } else {

    chain.l[i] <- lold

  }

  lcurrent <- chain.l[i]

  ## UPDATE sig2 ###############################################################
  sig2new <- rgamma(1, shape = sig2old/2, scale = 2)
  logHR <- logposterior(l = lcurrent, sig2 = sig2new,
                         ystar = ystarold, yprime = yprimeold,
                         given = given) -
    logposterior(l = lcurrent, sig2 = sig2old,
                  ystar = ystarold, yprime = yprimeold,
                  given = given)

  cat("sig2new logHR:", logHR, "\n")

  if ( !is.nan(logHR) & !is.na(logHR) ) {

    if (logHR > log(runif(1))) {
      chain.sig2[i] <- sig2new
    } else {
      chain.sig2[i] <- sig2old
    }

  } else {
    chain.sig2[i] <- sig2old
  }

  sig2current <- chain.sig2[i]

  ## UPDATE ystaryprime ########################################################
  S11 <- covMatrix(X = xstar, sig2 = sig2current,
                   covar.fun = r.matern, l = lcurrent)
  S22 <- covMatrix(X = xprime, sig2 = sig2current,
                   covar.fun = r.matern2, l = lcurrent)
  # CAREFUL SEE PAPER FOR ARG. ORDER
  S21 <- covMatrix(X = xprime, X2 = xstar, sig2 = sig2current,
                   covar.fun = r.matern1, l = lcurrent)
  S <- rbind(cbind(S11, t(S21)),
             cbind(S21, S22))

  #
  # Syy <- covMatrix(X = x, sig2 = sig2, covar.fun = r.matern, l = l)

  ystaryprimenew <- rmvnorm(1, mean = c(ystarold, yprimeold), sigma = v2 * S)
  ystarnew <- ystaryprimenew[1:length(xstar)]
  yprimenew <- ystaryprimenew[(length(xstar)+1):length(ystaryprimenew)]

  logHR <- logposterior(l = lcurrent, sig2 = sig2current,
                         ystar = ystarnew, yprime = yprimenew,
                         S = S,
                         given = given) -
    logposterior(l = lcurrent, sig2 = sig2current,
                  ystar = ystarold, yprime = yprimeold,
                  given = given)

  cat("ystarynew logHR:", logHR, "\n")

  if ( !is.nan(logHR) & !is.na(logHR) ) {

    if (logHR > log(runif(1))) {
      chain.ystar[,i] <- ystarnew
      chain.yprime[,i] <- yprimenew
    } else {
      chain.ystar[,i] <- ystarold
      chain.yprime[,i] <- yprimeold
    }

  } else {
    chain.ystar[,i] <- ystarold
    chain.yprime[,i] <- yprimeold
  }

  cat("Finished iteration:", i, "\n \n")

}
#)
# saveRDS(out, file = "/Users/ggh/Git/masters/eta0out.rds")

# after MH within Gibbs
chains <- list(eta.init = eta.init,
               chain.l = chain.l,
               chain.sig2 = chain.sig2,
               chain.ystar = chain.ystar,
               chain.yprime = chain.yprime)
saveRDS(chains, file = "/Users/ggh/Git/masters/eta0chains.rds")


chain.l
chain.sig2
cbind(eta.init$yprime, chain.yprime)

head(cbind(eta.init$ystar, chain.ystar))
pred.ystar <- apply(chain.ystar, MARGIN = 1, FUN = mean)
cbind(eta.init$ystar, pred.ystar)


