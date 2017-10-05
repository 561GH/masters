### Sampling functions #########################################################
library(foreach)
# return:
#     current: either accepted proposal, or stay at old for this iteration
#     accepted: 1 (accepted) or 0 (not accepted)

# ==============================================================================
update.l <- function(lold,
                     sig2, ystar, yprime,
                     given,
                     step) {

  # these stay the same if proposal not accepted
  lcurrent <- lold
  accepted <- 0

  lnew <- rmvn(1, m = lold, S = step)
  while (sum(lnew < 0) != 0) {  # ensures that all l aren't negative
    lnew <- rmvn(1, m = lold, S = step)
  }

  # TODO: check proposal with shirin
  # delta <- rnorm(1, 0, sd = step.l)
  # lnew <- lold + delta

  num <- logposterior2(l = lnew,
                       sig2 = sig2, ystar = ystar, yprime = yprime,
                       given = given) +
    log( dmvn(lold, m = lnew, S = step) )#/ pnorm(0.05, mean = lnew, sd = v1) )
  den <- logposterior2(l = lold,
                       sig2 = sig2, ystar = ystar, yprime = yprime,
                       given = given) +
    log( dmvn(lnew, m = lold, S = step) )#/ pnorm(0.05, mean = lold, sd = v1) )
  logHR <- num - den

  #cat("lnew logHR:", logHR)

  if ( !is.nan(logHR) & !is.na(logHR) ) {  # when have -Inf - Inf get NaN
    if ( logHR > log(runif(1)) )  {
      lcurrent <- lnew
      accepted <- 1
    }
  }

  return(list(lcurrent = lcurrent,
              accepted = accepted))

}

# ==============================================================================
update.sig2 <- function(l, sig2old,
                        ystar, yprime,
                        given) {

  # these stay the same if proposal not accepted
  sig2current <- sig2old
  accepted <- 0


  sig2new <- rchisq(1, df = sig2old)
  while (sig2new <= 0.5) {  # arbitrary minimum following shirin
    sig2new <- rchisq(1, df = sig2old)
  }

  num <- logposterior2(l = l, sig2 = sig2new,
                       ystar = ystar, yprime = yprime,
                       given = given) +
    log( pchisq(sig2new, df = sig2old) ) #/ pchisq(1.7, df = sig2old) )
  den <- logposterior2(l = l, sig2 = sig2old,
                       ystar = ystar, yprime = yprime,
                       given = given) +
    log( pchisq(sig2old, df = sig2new) ) #/ pchisq(1.7, df = sig2new) )
  logHR <- num - den

  #cat("sig2new logHR:", logHR)

  if ( !is.nan(logHR) & !is.na(logHR) ) {

    if (logHR > log(runif(1))) {
      sig2current <- sig2new
      accepted <- 1
    }

  }

  return(list(sig2current = sig2current,
              accepted = accepted))

}

# ==============================================================================
update.ystaryprime <- function(l, sig2, ystarold, yprimeold,
                               tau,
                               given,
                               step) {

  x <- given$x
  xstar <- given$xstar
  xprime <- given$xprime
  y <- given$y

  # these stay the same if proposal not accepted
  ystarcurrent <- ystarold
  yprimecurrent <- yprimeold
  accepted <- 0

  # did not need to do metropolis for eta_0 (particle initialisation)
  # because when tau0 = 0 (no constraint), we can sample from the
  # (ystar, ynew) | l, sig2 directly

  nugget <- 0
  pred.parameters <- predictiveMeanCov2(given = given,
                                        l = l, sig2 = sig2)
  m <- pred.parameters$mu.xstarprime
  tau2.xstarprime <- pred.parameters$cov.xstarprime +
    diag(rep(nugget, nrow(xstar) + nrow(xprime)))

  # trying to make an actual correlation matrix
  # V <- diag(diag(tau2.xstarprime))
  # Vsqrt.inv <- solve( (sqrt(V)) )
  # corr <- Vsqrt.inv %*% tau2.xstarprime %*% Vsqrt.inv

  # wtf don't get new samples
  # if (tau2.xstarprime[1,1] <= 1e-4) {
  #   ystarnew <- ystarold
  #   yprimenew <- yprimeold
  #
  #   return(list(ystarcurrent = ystarcurrent,
  #               yprimecurrent = yprimecurrent,
  #               accepted = NA,
  #               prob = NA,
  #               ystar2old = ystarold[2],
  #               ystar2new = ystarnew[2],
  #               yprime2old = yprimeold[2],
  #               yprime2new = yprimenew[2]))
  # }

  cat("Actually sampled wow.")

  # we have a symmetric proposal distribution
  # i.e. doesn't depend on previous ystar yprime value
  # ystaryprimenew <- rmvnorm(1, mean = c(ystarold, yprimeold),
  #                           sigma =  step * tau2.xstarprime) # qqq

  ystaryprimenew <- rmvn(1, m = c(ystarold, yprimeold),
                            S =  step * tau2.xstarprime)
  ystarnew <- ystaryprimenew[,1:nrow(xstar)]
  yprimenew <- ystaryprimenew[,-(1:nrow(xstar))]

  # constraint part is from likelihood
  num <- dmvn(c(ystarnew, yprimenew),
                 m = m, S = tau2.xstarprime, log = TRUE) +
    sum( log( pnorm(yprimenew * tau) ) )

  den <- dmvn(c(ystarold, yprimeold),
                 m = m, S = tau2.xstarprime, log = TRUE) +
    sum( log( pnorm(yprimeold * tau) ) )

  logHR <- num - den

  #cat("ystaryprime num, den, logHR:", num, den, logHR, "\n")
  #cat("ystaryprime logHR:", logHR)

  if ( !is.nan(logHR) & !is.na(logHR) ) {
    if ( logHR > log(runif(1)) ) {
        ystarcurrent <- ystarnew
        yprimecurrent <- yprimenew
        accepted <- 1
    }
  }

  return(list(ystarcurrent = ystarcurrent,
              yprimecurrent = yprimecurrent,
              accepted = accepted,
              prob = exp( logHR ),# qqq
              ystar2old = ystarold[2],
              ystar2new = ystarnew[2],
              yprime2old = yprimeold[2],
              yprime2new = yprimenew[2]))

}

# ==============================================================================

combine.particles <- function(out1, out2) {

  return( list(l = rbind(out1$l, out2$l),
               sig2 = c(out1$sig2, out2$sig2),
               ystar = rbind(out1$ystar, out2$ystar),
               yprime = rbind(out1$yprime, out2$yprime),
               acceptances = rbind(out1$acceptances, out2$acceptances)) )

}

update.all <- function(accepted = NULL,
                       steps,
                       particles.l.old,
                       particles.sig2.old,
                       particles.ystar.old,
                       particles.yprime.old,
                       tau) {

  # TODO: acceptance rates what
  accepted.l <- accepted.sig2 <- accepted.ystaryprime <- 0

  step.l <- steps$step.l
  step.ystaryprime <- steps$step.ystaryprime

  foreach(i = 1:N, .combine = "combine.particles") %dopar% {

    lold <- particles.l.old[i,]
    sig2old <- particles.sig2.old[i]
    ystarold <- particles.ystar.old[i,] # full ith particle for ystar as a vector
    yprimeold <- particles.yprime.old[i,] # full ith particle for yprime as a vector

    ## current: either accepted proposal, or stay at old for this iteration

    ### UPDATE l ###############################################################

    lupdate <- update.l(lold = lold,
                        sig2 = sig2old, ystar = ystarold, yprime = yprimeold,
                        given = given, step = step.l)

    lcurrent <- lupdate$lcurrent
    accepted.l <- accepted.l + lupdate$accepted

    #cat("\t \t \t -> acceptance rate for l:", accepted.l / i, "\n")

    ### UPDATE sig2 ############################################################

    sig2update <- update.sig2(l = lcurrent, sig2old = sig2old,
                              ystar = ystarold, yprime = yprimeold,
                              given = given)

    sig2current <- sig2update$sig2current
    accepted.sig2 <- accepted.sig2 + sig2update$accepted

    #cat("\t \t \t -> acceptance rate for sig2:", accepted.sig2 / i, "\n")

    ### UPDATE ystaryprime #####################################################
    # did not need to do metropolis for eta_0 (particle initialisation)
    # because when tau0 = 0 (no constraint), we can sample from the
    # (ystar, ynew) | l, sig2 directly

    ystaryprimeupdate <- update.ystaryprime(l = lcurrent, sig2 = sig2current,
                              ystar = ystarold, yprime = yprimeold,
                              given = given, step = step.ystaryprime,
                              tau = tau)

    ystarcurrent <- ystaryprimeupdate$ystarcurrent
    yprimecurrent <- ystaryprimeupdate$yprimecurrent
    accepted.ystaryprime <- accepted.ystaryprime + ystaryprimeupdate$accepted
    prob <- ystaryprimeupdate$prob  # qqq
    ystar2old <- ystaryprimeupdate$ystar2old
    ystar2new <- ystaryprimeupdate$ystar2new
    yprime2old <- ystaryprimeupdate$yprime2old
    yprime2new <- ystaryprimeupdate$yprime2new

    ### UPDATE chains ##########################################################

    list(l = lcurrent,
         sig2 = sig2current,
         ystar = ystarcurrent,
         yprime = yprimecurrent,
         acceptances =
           cbind(accepted.l = accepted.l,
                 accepted.sig2 = accepted.sig2,
                 accepted.ystaryprime = accepted.ystaryprime,
                 prob.ystaryprime = prob, # qqq
                 ystar2old = ystar2old,
                 ystar2new = ystar2new,
                 yprime2old = yprime2old,
                 yprime2new = yprime2new))

  }

}
