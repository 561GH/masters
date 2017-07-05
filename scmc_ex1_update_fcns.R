### Sampling functions #########################################################
library(foreach)
# return:
#     current: either accepted proposal, or stay at old for this iteration
#     accepted: 1 (accepted) or 0 (not accepted)

# ==============================================================================
update.l <- function(lold,
                     sig2, ystar, yprime,
                     given = given,
                     step = step.l) {

  # these stay the same if proposal not accepted
  lcurrent <- lold
  accepted <- 0

  lnew <- rnorm(1, mean = lold, sd = step.l)  # TODO: check proposal with shirin

  if (lnew > 0) {

    logHR <- logposterior(l = lnew,
                          sig2 = sig2, ystar = ystar, yprime = yprime,
                          given = given) -
      logposterior(l = lold,
                   sig2 = sig2, ystar = ystar, yprime = yprime,
                   given = given)

    #cat("lnew logHR:", logHR)

    if ( !is.nan(logHR) & !is.na(logHR) ) {  # when have -Inf - Inf get NaN
      if ( logHR > log(runif(1)) )  {
        lcurrent <- lnew
        accepted <- 1
      }
    }

  }

  return(list(lcurrent = lcurrent,
              accepted = accepted))

}

# ==============================================================================
update.sig2 <- function(l, sig2old,
                        ystar, yprime,
                        given = given,
                        step = step.sig2) {

  # these stay the same if proposal not accepted
  sig2current <- sig2old
  accepted <- 0

  sig2new <- rchisq(1, df = sig2old) #rgamma(1, shape = sig2old/2, scale = 2)

  logHR <- logposterior(l = l, sig2 = sig2new,
                        ystar = ystar, yprime = yprime,
                        given = given) -
    logposterior(l = l, sig2 = sig2old,
                 ystar = ystar, yprime = yprime,
                 given = given)

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
                        given = given,
                        step = step.ystarynew) {

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
  nugget <- 1e-6  # -7 caused 39/100 warnings

  ## fix this solve; check the nugget;
  # OR check with another inversion method
  Rinv <- solve(covMatrix(X = x, sig2 = sig2,
                          covar.fun = r.matern, l = l) +
                  diag(nugget, nrow(x)))

  S11 <- covMatrix(X = xstar, sig2 = sig2, covar.fun = r.matern, l = l)
  S22 <- covMatrix(X = xprime, sig2 = sig2, covar.fun = r.matern2, l = l)
  S21 <- covMatrix(X = xprime, X2 = xstar,  # CAREFUL SEE PAPER FOR ARG. ORDER
                   sig2 = sig2, covar.fun = r.matern1, l = l)
  R.xstarxprime <- rbind(cbind(S11, t(S21)),
                         cbind(S21, S22)) +
    diag(nugget, nrow(xstar) + nrow(xprime))

  # CAREFUL SEE PAPER FOR ARG. ORDER
  r.xstarprime.x <- rbind(covMatrix(X = xstar, X2 = x, sig2 = sig2,
                                    covar.fun = r.matern, l = l),
                          covMatrix(X = xprime, X2 = x, sig2 = sig2,
                                    covar.fun = r.matern1, l = l))
  mu.xstarprime <- r.xstarprime.x %*% Rinv %*% y

  # because covMatrix has sig2 in it, have an extra one multiplying in???
  tau2.xstarprime <- R.xstarxprime -
    r.xstarprime.x %*% Rinv %*% t(r.xstarprime.x) #/ sig2current #???

  # MAKE SURE tau2.xstarprime IS SYMMETRIC
  tau2.xstarprime <- ( tau2.xstarprime + t(tau2.xstarprime) ) / 2

  # MAKE SURE tau2.xstarprime IS POSTIVE SEMI-DEFINITE
  # eigen(tau2.xstarprime)$values  # changed nugget to 10e-6 from -8

  ystaryprimenew <- rmvnorm(1, mean = mu.xstarprime, sigma = tau2.xstarprime)
  ystarnew <- ystaryprimenew[,1:nrow(xstar)]
  yprimenew <- ystaryprimenew[,-(1:nrow(xstar))]

  logHR <- logposterior(l = l, sig2 = sig2, ystar = ystarnew, yprime = yprimenew,
                        given = given) -
    logposterior(l = l, sig2 = sig2, ystar = ystarold, yprimeold,
                 given = given)

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
              accepted = accepted))

}

# ==============================================================================

combine.particles <- function(out1, out2) {

  return( list(l = c(out1$l, out2$l),
               sig2 = c(out1$sig2, out2$sig2),
               ystar = rbind(out1$ystar, out2$ystar),
               yprime = rbind(out1$yprime, out2$yprime),
               acceptances = rbind(out1$acceptances, out2$acceptances)) )

}

update.all <- function(accepted = NULL,
                       steps = NULL,
                       particles.l.old,
                       particles.sig2.old,
                       particles.ystar.old,
                       particles.yprime.old) {

  # TODO: acceptance rates what
  accepted.l <- accepted.sig2 <- accepted.ystaryprime <- 0

  foreach(i = 1:N, .combine = "combine.particles") %dopar% {

    lold <- particles.l.old[i]
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
                              given = given, step = step.sig2)

    sig2current <- sig2update$sig2current
    accepted.sig2 <- accepted.sig2 + sig2update$accepted

    #cat("\t \t \t -> acceptance rate for sig2:", accepted.sig2 / i, "\n")

    ### UPDATE ystaryprime #####################################################
    # did not need to do metropolis for eta_0 (particle initialisation)
    # because when tau0 = 0 (no constraint), we can sample from the
    # (ystar, ynew) | l, sig2 directly

    ystaryprimeupdate <- update.ystaryprime(l = lcurrent, sig2 = sig2current,
                              ystar = ystarold, yprime = yprimeold,
                              given = given, step = step.ystaryprime)

    ystarcurrent <- ystaryprimeupdate$ystarcurrent
    yprimecurrent <- ystaryprimeupdate$yprimecurrent
    accepted.ystaryprime <- accepted.ystaryprime + ystaryprimeupdate$accepted

    ### UPDATE chains ##########################################################

    list(l = lcurrent,
         sig2 = sig2current,
         ystar = ystarcurrent,
         yprime = yprimecurrent,
         acceptances =
           cbind(accepted.l = accepted.l,
                 accepted.sig2 = accepted.sig2,
                 accepted.ystaryprime = accepted.ystaryprime))

  }

}
