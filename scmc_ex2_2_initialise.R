# SCMC initialise particles

## run setup ###################################################################

# for stepping through manually
#N <- 100
#eta.init <- initial.values
#v1 <- diag(c(0.07, 0.07))
#i <- 1

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

  chain.l <- matrix(NA, nrow = N, ncol = ncol(x))
  chain.sig2 <- rep(NA, N)
  chain.ystar <- matrix(NA, nrow = N, ncol = nrow(xstar))
  chain.yprime <- matrix(NA, nrow = N, ncol = nrow(xprime))

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
      lold <- chain.l[i-1,]
      sig2old <- chain.sig2[i-1]
      ystarold <- chain.ystar[i-1,]
      yprimeold <- chain.yprime[i-1,]
    }

    chain.l[i,] <- lold
    chain.sig2[i] <- sig2old
    chain.ystar[i,] <- ystarold
    chain.yprime[i,] <- yprimeold

    ## UPDATE l ################################################################
    lnew <- rmvn(1, m = lold, S = v1)
    while (sum(lnew < 0.05) != 0) {  # ensures that all l aren't too small
      lnew <- rmvn(1, m = lold, S = v1)
    }

    num <- logposterior2(l = lnew,
                        sig2 = sig2old, ystar = ystarold, yprime = yprimeold,
                        given = given) +
      log( dmvn(lold, m = lnew, S = v1) )#/ pnorm(0.05, mean = lnew, sd = v1) )
    den <- logposterior2(l = lold,
                        sig2 = sig2old, ystar = ystarold, yprime = yprimeold,
                        given = given) +
      log( dmvn(lnew, m = lold, S = v1) )#/ pnorm(0.05, mean = lold, sd = v1) )
    logHR <- num - den

    #cat("lnew logHR:", logHR)

    if ( !is.nan(logHR) & !is.na(logHR) ) {  # when have -Inf - Inf get NaN
      if ( logHR > log(runif(1)) )  {
        chain.l[i,] <- lnew
        accepted.l[i] <- 1
      }
    }

    lcurrent <- chain.l[i,]
    #cat("\t \t \t -> acceptance rate for l:", sum(accepted.l) / i, "\n")

    ## UPDATE sig2 ###############################################################
    sig2new <- rchisq(1, df = sig2old)
    while (sig2new <= 1.7) {  # arbitrary minimum following shirin
      sig2new <- rchisq(1, df = sig2old)
    }

    num <- logposterior2(l = lcurrent, sig2 = sig2new,
                        ystar = ystarold, yprime = yprimeold,
                        given = given) +
      log( pchisq(sig2new, df = sig2old) ) #/ pchisq(1.7, df = sig2old) )
    den <- logposterior2(l = lcurrent, sig2 = sig2old,
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
    nugget <- 0
    pred.parameters <- predictiveMeanCov2(given = given,
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
# burn <- 500
# N <- 1000
# init <- eta0(eta.init = initial.values,
#              given = given,  # data, locations, obs, etc.)
#              N = burn + N, # particles
#              v1 = diag(c(0.1, 0.1))) # step size for l proposal
#
# particles.l <- init$chain.l[-(1:burn),]
# particles.sig2 <- init$chain.sig2[-(1:burn)]
#
# hist(particles.l[,1])
# hist(particles.l[,2])
# hist(particles.sig2)
#
# plot(particles.l[,1], type = "l")
# plot(particles.l[,2], type = "l")
# plot(particles.sig2, type = "l")
#
# # compare with true values
# particles.yprime <- init$chain.yprime[-(1:burn),]
# tail(particles.yprime)
# true.values$yprime
#
# # par(mfrow = c(floor(sqrt(ncol(particles.yprime))),
# #               ceiling(sqrt(ncol(particles.yprime)))))
# for (j in 1:ncol(particles.yprime)) {
#   pdf(file = paste("yprime", j, ".pdf", sep = ""))
#   hist(particles.yprime[,j])
#   abline(v = true.values$yprime[j], col = "red")
#   dev.off()
# }
# # par(mfrow = c(1, 1))
#
# ###
# particles.ystar <- init$chain.ystar[-(1:burn),]
# tail(particles.ystar)
# true.values$ystar
# par(mfrow = c(floor(sqrt(ncol(particles.ystar))),
#               ceiling(sqrt(ncol(particles.ystar)))))
# for (j in 1:ncol(particles.ystar)) {
#   hist(particles.ystar[,j])
#   abline(v = true.values$ystar[j], col = "red")
# }
# par(mfrow = c(1, 1))
