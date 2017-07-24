# SCMC initialise particles

## run setup ###################################################################

## start MH within gibbs #######################################################
# i.e. unconstrained SCMC

eta0 <- function(eta.init,  # initial values for l, sig2, ystar, yprime is list
                 given,  # data, locations, obs, etc.
                 N, # iterations
                 v1) {  # step size for l proposal

  x <- given$x
  xstar <- given$xstar
  xprime <- given$xprime
  y <- given$y

  chain.l <- rep(NA, N)
  chain.sig2 <- rep(NA, N)
  chain.ystar <- matrix(NA, nrow = N, ncol = length(xstar))
  chain.yprime <- matrix(NA, nrow = N, ncol = length(xprime))

  ### Track acceptance rate
  accepted.l <- accepted.sig2 <- rep(0, N)

  for (i in 1:N) {

    cat("Starting iteration:", i, "out of", N, "\n")

    if (i == 1) {
      lold <- eta.init$l
      sig2old <- eta.init$sig2
      ystarold <- eta.init$ystar
      yprimeold <- eta.init$yprime
    } else {
      lold <- chain.l[i-1]
      sig2old <- chain.sig2[i-1]
      ystarold <- chain.ystar[i-1,]
      yprimeold <- chain.yprime[i-1,]
    }

    chain.l[i] <- lold
    chain.sig2[i] <- sig2old
    chain.ystar[i,] <- ystarold
    chain.yprime[i,] <- yprimeold

    ## UPDATE l ##################################################################
    lnew <- rnorm(1, mean = lold, sd = v1)
    while (lnew < 0.05) {  # ensures that l isn't too stupidly small
      lnew <- rnorm(1, mean = lold, sd = v1)
    }

    num <- logposterior(l = lnew,
                        sig2 = sig2old, ystar = ystarold, yprime = yprimeold,
                        given = given) +
      log( pnorm(lold, mean = lnew, sd = v1) )#/ pnorm(0.05, mean = lnew, sd = v1) )
    den <- logposterior(l = lold,
                        sig2 = sig2old, ystar = ystarold, yprime = yprimeold,
                        given = given) +
      log( pnorm(lnew, mean = lold, sd = v1) )#/ pnorm(0.05, mean = lold, sd = v1) )
    logHR <- num - den

    #cat("lnew logHR:", logHR)

    if ( !is.nan(logHR) & !is.na(logHR) ) {  # when have -Inf - Inf get NaN
      if ( logHR > log(runif(1)) )  {
        chain.l[i] <- lnew
        accepted.l[i] <- 1
      }
    }

    lcurrent <- chain.l[i]
    #cat("\t \t \t -> acceptance rate for l:", sum(accepted.l) / i, "\n")

    ## UPDATE sig2 ###############################################################
    sig2new <- rchisq(1, df = sig2old)
    while (sig2new <= 1.7) {  # arbitrary minimum following shirin
      sig2new <- rchisq(1, df = sig2old)
    }

    num <- logposterior(l = lcurrent, sig2 = sig2new,
                        ystar = ystarold, yprime = yprimeold,
                        given = given) +
      log( pchisq(sig2new, df = sig2old) ) #/ pchisq(1.7, df = sig2old) )
    den <- logposterior(l = lcurrent, sig2 = sig2old,
                        ystar = ystarold, yprime = yprimeold,
                        given = given) +
      log( pchisq(sig2old, df = sig2new) ) #/ pchisq(1.7, df = sig2new) )
    logHR <- num - den

    #cat("sig2new logHR:", logHR)

    if ( !is.nan(logHR) & !is.na(logHR) ) {

      if (logHR > log(runif(1))) {
        chain.sig2[i] <- sig2new
        accepted.sig2[i] <- 1
      }

    }

    sig2current <- chain.sig2[i]
    #cat("\t \t \t -> acceptance rate for sig2:", sum(accepted.sig2) / i, "\n")

    ## UPDATE ystaryprime ########################################################
    nugget <- 1e-6

    pred.parameters <- predictiveMeanCov(given = given,
                                         l = lcurrent, sig2 = sig2current)
    mu.xstarprime <- pred.parameters$mu.xstarprime
    tau2.xstarprime <- pred.parameters$cov.xstarprime +
      diag(rep(nugget, nrow(given$xstar) + nrow(given$xprime)))

    # MAKE SURE tau2.xstarprime IS POSTIVE SEMI-DEFINITE
    # eigen(tau2.xstarprime)$values

    ystaryprimenew <- rmvn(1, m = mu.xstarprime, S = tau2.xstarprime)
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
# burn <- 1
# N <- 100
# init <- eta0(eta.init = list(l = 0.5, sig2 = 5,
#                              ystar = ytrue(given$xstar) - mean(ytrue(given$x)),
#                              yprime = 20 / (20 * given$xprime) ),
#              given = given,  # data, locations, obs, etc.)
#              N = burn + N, # particles
#              v1 = 0.07) # step size for l proposal
#
# hist(init$chain.l[-(1:burn)])
# hist(init$chain.sig2[-(1:burn)])
#
#
# particles.yprime <- init$chain.yprime[-(1:burn),]
# plot(x = given$xprime, y = 20 / (20 * given$xprime), pch = 16, cex = 2)
# for (n in (N-50):N) {
#   points(x = given$xprime, y = particles.yprime[n,], type = "l", col = n)
# }
# points(x = given$xprime, y = apply(particles.yprime, 2, mean), pch = 16, col = "red")
# tail(particles.yprime)
#
# particles.ystar <- init$chain.ystar[-(1:burn),]
# tail(particles.ystar)
# plot(x = given$x, y = given$y, pch = 16, cex = 2)
# for (n in (N-50):N) {
#   points(x = given$xstar, y = particles.ystar[n,], type = "l", col = "red")
# }
#
# particles.l <- init$chain.l[-(1:burn)]
# particles.sig2 <- init$chain.sig2[-(1:burn)]
#
# plot(particles.l, type = "l")
# plot(particles.sig2, type = "l")
#
# ###
# par(mfrow = c(2, 1))
# plot(x = 1:N, y = cumsum(init$accepted.l[-(1:burn)]) / 1:N, pch = 16,
#      ylab = "rate l")
# abline(h = 0.2); abline(h = 0.1)
# plot(x = 1:N, y = cumsum(init$accepted.sig2[-(1:burn)]) / 1:N, pch = 16,
#      ylab = "rate sig2")
# abline(h = 0.2); abline(h = 0.1)
# par(mfrow = c(1, 1))
