# SCMC initialise particles

## run setup ###################################################################
library(mvtnorm)
ytrue <- function(x) {log(20 * x + 1)}  # x > -1/20; monotone increasing function

# x = inputs where fcn evaluated
# xstar = prediction set inputs
# xprime = derivative set inputs
# true/known model values
given <- list(x = cbind(c(0, 0.1, 0.2, 0.3, 0.4, 0.9, 1)),
              xstar = cbind(seq(0, 1, length.out = 50)),
              xprime = cbind(c(0.42, 0.47, 0.52, 0.57, 0.62,
                               0.67, 0.72, 0.77, 0.82, 0.87)),
              y = ytrue(c(0, 0.1, 0.2, 0.3, 0.4, 0.9, 1)) -
                mean(ytrue(c(0, 0.1, 0.2, 0.3, 0.4, 0.9, 1))) )

## start MH within gibbs #######################################################
# i.e. unconstrained SCMC

eta0 <- function(eta.init,  # initial values for l, sig2, ystar, yprime is list
                 given,  # data, locations, obs, etc.
                 N = 100, # iterations
                 v1 = 0.01,  # step size for l proposal
                 v2 = 3) {# step size for sig2 proposal

  x <- given$x
  xstar <- given$xstar
  xprime <- given$xprime
  y <- given$y

  chain.l <- rep(NA, N)
  chain.sig2 <- rep(NA, N)
  chain.ystar <- matrix(NA, nrow = N, ncol = length(xstar))
  chain.yprime <- matrix(NA, nrow = N, ncol = length(xprime))

  ### Track acceptance rate
  accepted.l <- accepted.sig2 <- 0

  for (i in 1:N) {

    cat("Starting iteration:", i, "out of", N, "\n")

    if (i == 1) {
      lold <- eta.init$l
      sig2old <- eta.init$sig2
      ystarold <- eta.init$ystar
      yprimeold <- eta.init$yprime
    } else {
      lold <- chain.l[i-1]
      sig2old <- chain.l[i-1]
      ystarold <- chain.ystar[i-1,]
      yprimeold <- chain.yprime[i-1,]
    }

    chain.l[i] <- lold
    chain.sig2[i] <- sig2old
    chain.ystar[i,] <- ystarold
    chain.yprime[i,] <- yprimeold

    ## UPDATE l ##################################################################
    lnew <- rnorm(1, mean = lold, sd = v1)

    if (lnew > 0) {

      logHR <- logposterior(l = lnew,
                            sig2 = sig2old, ystar = ystarold, yprime = yprimeold,
                            given = given) -
        logposterior(l = lold,
                     sig2 = sig2old, ystar = ystarold, yprime = yprimeold,
                     given = given)

      #cat("lnew logHR:", logHR)

      if ( !is.nan(logHR) & !is.na(logHR) ) {  # when have -Inf - Inf get NaN
        if ( logHR > log(runif(1)) )  {
          chain.l[i] <- lnew
          accepted.l <- accepted.l + 1
        }
      }

    }

    lcurrent <- chain.l[i]
    cat("\t \t \t -> acceptance rate for l:", accepted.l / i, "\n")

    ## UPDATE sig2 ###############################################################
    sig2new <- rchisq(1, df = sig2old) #rgamma(1, shape = sig2old/2, scale = 2)
    logHR <- logposterior(l = lcurrent, sig2 = sig2new,
                          ystar = ystarold, yprime = yprimeold,
                          given = given) -
      logposterior(l = lcurrent, sig2 = sig2old,
                   ystar = ystarold, yprime = yprimeold,
                   given = given)

    #cat("sig2new logHR:", logHR)

    if ( !is.nan(logHR) & !is.na(logHR) ) {

      if (logHR > log(runif(1))) {
        chain.sig2[i] <- sig2new
        accepted.sig2 <- accepted.sig2 + 1
      }

    }

    sig2current <- chain.sig2[i]
    cat("\t \t \t -> acceptance rate for sig2:", accepted.sig2 / i, "\n")

    ## UPDATE ystaryprime ########################################################
    nugget <- 1e-6
    R <- covMatrix(X = x, sig2 = sig2current,
                   covar.fun = r.matern, l = lcurrent) +
      diag(nugget, nrow(x))
    # Linv <- solve( t(chol(R)) )  # recall need transpose to match std chol def
    # Rinv <- t(Linv) %*% Linv
    Rinv <- solve(R)

    S11 <- covMatrix(X = xstar, sig2 = sig2current, covar.fun = r.matern, l = lcurrent)
    S22 <- covMatrix(X = xprime, sig2 = sig2current, covar.fun = r.matern2, l = lcurrent)
    S21 <- covMatrix(X = xprime, X2 = xstar,  # CAREFUL SEE PAPER FOR ARG. ORDER
                     sig2 = sig2current, covar.fun = r.matern1, l = lcurrent)
    R.xstarxprime <- rbind(cbind(S11, t(S21)),
                           cbind(S21, S22)) +
      diag(nugget, nrow(xstar) + nrow(xprime))

    # CAREFUL SEE PAPER FOR ARG. ORDER
    r.xstarprime.x <- rbind(covMatrix(X = xstar, X2 = x, sig2 = sig2current,
                                      covar.fun = r.matern, l = lcurrent),
                            covMatrix(X = xprime, X2 = x, sig2 = sig2current,
                                      covar.fun = r.matern1, l = lcurrent))
    mu.xstarprime <- r.xstarprime.x %*% Rinv %*% y

    # because covMatrix has sig2 in it, have an extra one multiplying in???
    tau2.xstarprime <- R.xstarxprime -
      r.xstarprime.x %*% Rinv %*% t(r.xstarprime.x) #/ sig2current #???

    # MAKE SURE tau2.xstarprime IS SYMMETRIC
    tau2.xstarprime <- ( tau2.xstarprime + t(tau2.xstarprime) ) / 2

    # MAKE SURE tau2.xstarprime IS POSTIVE SEMI-DEFINITE
    # eigen(tau2.xstarprime)$values

    ystaryprimenew <- rmvnorm(1, mean = mu.xstarprime, sigma = tau2.xstarprime)
    ystarnew <- ystaryprimenew[,1:nrow(xstar)]
    yprimenew <- ystaryprimenew[,-(1:nrow(xstar))]

    chain.ystar[i,] <- ystarnew
    chain.yprime[i,] <- yprimenew

    #cat("Finished iteration:", i, "\n \n")

  }

  # after MH within Gibbs
  chains <- list(eta.init = eta.init,
                 chain.l = chain.l,
                 chain.sig2 = chain.sig2,
                 chain.ystar = chain.ystar,
                 chain.yprime = chain.yprime,
                 accepted.l = accepted.l,
                 accepted.sig2 = accepted.sig2)
  #saveRDS(chains, file = "/Users/ggh/Git/masters/eta0chains.rds")

  return(chains)
}

################################################################################
# # test
# burn <- 7
# N <- 5
# init <- eta0(eta.init = list(l = 5.7, sig2 = 1,
#                              ystar = ytrue(given$xstar) - mean(ytrue(given$x)),
#                              yprime = 20 / (20 * given$xprime) ),
#              given = given,  # data, locations, obs, etc.)
#              N = burn + N, # particles
#              v1 = 0.01, # step size for l proposal
#              v2 = 3) # step size for sig2 proposal
