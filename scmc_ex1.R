# SCMC toy example

## setup
ytrue <- function(x) {log(20 * x + 1)}  # x > -1/20; monotone increasing function

# x = inputs where fcn evaluated
# xstar = prediction set inputs
# xprime = derivative set inputs
# true/known model values
given <- list(x = cbind(c(0, 0.1, 0.2, 0.3, 0.4, 0.9, 1)),
              xstar = cbind(seq(0, 1, length.out = 50)),
              xprime = cbind(c(0.42, 0.47, 0.52, 0.57, 0.62,
                               0.67, 0.72, 0.77, 0.82, 0.87)),
              y = ytrue(c(0, 0.1, 0.2, 0.3, 0.4, 0.9, 1)) )

TT <- 20  # total time
tauTT <- 10e-6
taus <- seq(0, tauTT, length.out = TT+1)  # TODO: sequence of taus???
taus <- taus[2:length(taus)]

N <- 3  # particles

################################################################################
# 1. initialise (ith row = ith particle) #######################################
################################################################################
# jth colunn = jth time step
# - store multidimensional things in a list (ystar, yprime)
particles.l <- matrix(NA, ncol = TT, nrow = N)
particles.sig2 <- matrix(NA, ncol = TT, nrow = N)
particles.ystar <- lapply(rep(NA, nrow(given$xstar)),
                          function(x) matrix(NA, ncol = 1, nrow = N) )
particles.yprime <- lapply(rep(NA, nrow(given$xprime)),
                           function(x) matrix(NA, ncol = 1, nrow = N) )

# TODO: initialise t0 i.e. tau0 particles from eta0
particles.l.t0 <- NULL
particles.sig2.t0 <- NULL
particles.ystar.t0 <- NULL
particles.yprime.t0 <- NULL

################################################################################
# 2. weights at time t = 1 #####################################################
################################################################################
## matrix of weights: jth column = weights for t-th iteration
W <- matrix(NA, ncol = TT, nrow = N)
W[,1] <- 1/N

################################################################################
# 3. looping through t = 1, ..., TT-1 ##########################################
################################################################################
step.l <- 0.01  # step size for l
step.sig2 <- 6  # step size for sig2
step.ystarynew <- 1  # tune this

# should let step size vary to provide decent acceptance rates (during burn in?)

### Sampling functions #########################################################
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

  lnew <- rnorm(1, mean = 5.7, sd = step.l)

  if (lnew > 0) {

    logHR <- logposterior(l = lnew,
                          sig2 = sig2, ystar = ystar, yprime = yprime,
                          given = given) -
      logposterior(l = lold,
                   sig2 = sig2, ystar = ystar, yprime = yprime,
                   given = given)

    cat("lnew logHR:", logHR)

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

  cat("sig2new logHR:", logHR)

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
  ystarnew <- ystaryprimenew[1:nrow(xstar),]
  yprimenew <- ystaryprimenew[-(1:nrow(xstar)),]

  logHR <- logposterior(l = l, sig2 = sig2, ystar = ystarnew, yprime = yprimenew,
                        given = given) -
    logposterior(l = l, sig2 = sig2, ystar = ystarold, yprimeold,
                 given = given)

  cat("ystaryprime logHR:", logHR)

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

update.all <- function(t = t, accepted,  # TODO: acceptance rates???
                       particles.l.t = particles.l.t,
                       particles.sig2.t = particles.sig2.t,
                       particles.ystar.t = particles.ystar.t,
                       particles.yprime.t = particles.yprime.t) { # TODO: fill this in

  foreach(i = 1:N, .combine = "combine.particles") %dopar% {

    lold <- particles.l.t[i]
    sig2old <- particles.sig2.t[i]
    # all the "parts" for a single particle are spread accross elements
    # in a list => need to get the ith entry from each element of the list
    # to get the full ith particle for ystar and yprime
    ystarold <- lapply(1:nrow(given$xstar),
                       function(x) (particles.ystar.t[[x]])[i])
    ystarold <- unlist(ystarold)  # full ith particle for ystar as a vector
    yprimeold <- lapply(1:nrow(given$xprime),
                       function(x) (particles.xprime.t[[x]])[i])
    yprimeold <- unlist(yprimeold)  # full ith particle for yprime as a vector

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

    ystaryprimeupdate <- update.ystaryprime(l = lcurrent, sig2old = sig2current,
                              ystar = ystarold, yprime = yprimeold,
                              given = given, step = step.ystaryprimenew)

    ystarcurrent <- ystaryprimeupdate$ystarcurrent
    yprimecurrent <- ystaryprimeupdate$yprimecurrent
    accepted.ystaryprime <- accepted.ystaryprime + ystaryprimeupdate$accepted

    ### UPDATE chains ##########################################################

    return( list(l = lcurrent,
                 sig2 = sig2current,
                 ystar = ystarcurrent,
                 yprime = yprimecurrent,
                 acceptances =
                    cbind(accepted.l = accepted.l,
                          accepted.sig2 = accepted.sig2,
                          accepted.ystaryprime = accepted.ystaryprime)) )

  }

}

for (t in 1:TT) {

  ## WEIGHT CALCULATION ========================================================
  ## AND t-1 (old) PARTICLES ===================================================

  if (t == 1) {

    # TODO: check if you really set particles_1 to particles_0
    particles.l.old <- particles.l.t <- particles.l.t0
    particles.sig2.old <- particles.sig2.t <- particles.sig2.t0
    particles.ystar.old <- particles.ystar.t <- particles.ystar.t0
    particles.yprime.old <- particles.yprime.t <- particles.yprime.t0

    # ==========================================================================
    yprimet <- apply(particles.yprime.old, MARGIN = 2, FUN = mean)

    wtildes <- prod( pnorm(tau[1] * yprimet) ) /
      prod( pnorm(0 * yprimet) )  # tau0 = 0

    W[,1] <- wtildes

  } else {

    # t - 1 ITERATION ==========================================================
    particles.l.old <- particles.l[,t-1]
    particles.sig2.old <- particles.sig2[,t-1]
    particles.ystar.old <- lapply(1:nrow(given$xstar),
                                function(x) (particles.ystar[[x]])[,t-1])
    particles.yprime.old <- lapply(1:nrow(given$xprime),
                                function(x) (particles.yprime[[x]])[,t-1])

    # t ITERATION ==============================================================
    particles.l.t <- particles.l[,t]
    particles.sig2.t <- particles.sig2[,t]
    # each list element is for a single input
    particles.ystar.t <- lapply(1:nrow(given$xstar),
                                function(x) (particles.ystar[[x]])[,t])
    particles.yprime.t <- lapply(1:nrow(given$xprime),
                                function(x) (particles.yprime[[x]])[,t])

    # ==========================================================================
    yprimet <- apply(particles.yprime.old, MARGIN = 2, FUN = mean)

    wtildes <- prod( pnorm(tau[t] * particles.yprime.old) ) /
      prod( pnorm(tau[t-1] * particles.yprime.old) )

    W[,t] <- W[,t-1] * wtildes
  }

  W[,t] <- W[,t] / sqrt( sum(W[,t]^2) )  # normalise weights

  ## RESAMPLING ================================================================
  ### effective sample size
  ESSt <- 1 / sum( (W[,t])^2 )

  if (ESSt <= N/2) {

    resample <- sample(1:N, replace = TRUE, prob = W[,t])

    particles.l.t <- particles.l.t[resample]
    particles.sig2.t <- particles.sig2.t[resample]
    particles.ystar.t <- lapply(1:nrow(given$xstar),
                                function(x) (particles.ystar.t[[x]])[resample])
    particles.yprime.t <- lapply(1:nrow(given$xprime),
                                function(x) (particles.yprime.t[[x]])[resample])

    W[,t] <- 1/N

  }

  ## SAMPLING ==================================================================

  ### t + 1 iteration preparation ==============================================
  particles.l.t <- particles.l
  particles.sig2.t <- particles.sig2
  particles.ystar.t <- particles.ystar
  particles.yprime.t <- particles.yprime

  ### Track acceptance rate ====================================================
  accepted.l <- accepted.sig2 <- accepted.ystaryprime <- 0

  ### Sample and update chains =================================================
  particles.new <- update.all()  # TODO: fill this in after update.all done (acceptance rates ???)
  acceptances <- particles.new$acceptances  # N by 3 matrix of acceptances

  particles.l[,t+1] <- particles.new$sig2
  particles.sig2[,t+1] <- particles.new$sig2current

  ystar.new <- particles.new$ystar
  yprime.new <- particles.new$yprime

  particles.ystar <- lapply(1:length(given$xstar),
                            function(x) cbind(particles.ystar[[x]], ystar.new[,x]))

  particles.yprime <- lapply(1:length(given$xprime),
                             function(x) cbind(particles.yprime[[x]], yprime.new[,x]))
  particles.ystarcurrent  # TODO: adding new particles to chain
                              #  particles.sig2 <- matrix(NA, ncol = TT, nrow = N)
                              #  particles.ystar <- lapply(rep(NA, nrow(given$xstar)),
                              #                            function(x) matrix(NA, ncol = TT, nrow = N) )
                              #  particles.yprime <- lapply(rep(NA, nrow(given$xprime)),
                              #                             function(x) matrix(NA, ncol = TT, nrow = N) )



}
