# SCMC first actual attempt

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

TT <- 20  # total time
tauTT <- 10e-6
taus <- seq(0, tauTT, length.out = TT+1)  # qqq for now
taus <- taus[2:length(taus)]

N <- 100  # particles

################################################################################
# 1. initialise (ith row = ith particle) #######################################
################################################################################
particles.l <- matrix(1, ncol = 1, nrow = N)
particles.sig2 <- matrix(1, ncol = 1, nrow = N)
particles.ystar <- matrix(NA, ncol = length(xstar), nrow = N)
particles.yprime <- matrix(NA, ncol = length(xprime), nrow = N)

################################################################################
# 2. weights at time t = 1 #####################################################
################################################################################
## matrix of weights: jth column = weights for t-th iteration
W <- matrix(NA, ncol = TT, nrow = N)
W[,1] <- 1/N

################################################################################
# 3. looping through t = 1, ..., TT-1 ##########################################
################################################################################
v1 <- 0.01  # step size for l
v2 <- 6  # step size for sig2

# should let v1 and v2 vary to provide decent acceptance rates

for (t in 1:TT) {

  ## WEIGHT CALCULATION ========================================================
  yprimet <- apply(particles.yprime, MARGIN = 2, FUN = mean)

  if (t == 1) {
    wtildes <- prod( pnorm(tau[1] * yprimet) ) /
      prod( pnorm(0 * yprimet) )  # tau0 = 0

    W[,1] <- wtildes

  } else {
    wtildes <- prod( pnorm(tau[t] * particles.yprime) ) /
      prod( pnorm(tau[t-1] * particles.yprime) )

    W[,t] <- W[,t-1] * wtildes
  }

  W[,t] <- W[,t] / sqrt( sum(W[,t]^2) )  # normalise weights

  ## RESAMPLING ================================================================
  ### effective sample size
  ESSt <- 1 / sum(W[,t])

  if (ESSt <= N/2) {

    resample <- sample(1:N, replace = TRUE, prob = W[,t])

    particles.l <- particles.l[resample]
    particles.sig2 <- particles.sig2[resample]
    particles.ystar <- particles.ystar[resample,]
    particles.yprime <- particles.yprime[resample,]

    W[,t] <- 1/N

  }

  ## SAMPLING ==================================================================

  particles.l.t <- particles.l
  particles.sig2.t <- particles.sig2
  particles.ystar.t <- particles.ystar
  particles.yprime.t <- particles.yprime

  ### Track acceptance rate
  accepted.l <- accepted.sig2 <- 0

  particles.t <-
    foreach(i = 1:N) %dopar% {  # .combine??

      lold <- particles.l[i]
      sig2old <- particles.sig2[i]
      ystarold <- particles.ystar[i,]
      yprimeold <- particles.yprime[i,]

      ### UPDATE l #############################################################
      lnew <- rnorm(1, mean = 5.7, sd = v1)

      if (lnew > 0) {

        logHR <- logposterior(l = lnew,
                              sig2 = sig2old, ystar = ystarold, yprime = yprimeold,
                              given = given) -
          logposterior(l = lold,
                       sig2 = sig2old, ystar = ystarold, yprime = yprimeold,
                       given = given)

        cat("lnew logHR:", logHR)

        if ( !is.nan(logHR) & !is.na(logHR) ) {  # when have -Inf - Inf get NaN
          if ( logHR > log(runif(1)) )  {
            particles.l.t[i] <- lnew
            accepted.l <- accepted.l + 1
          }
        }

      }

      lcurrent <- particles.l.t[i]
      #cat("\t \t \t -> acceptance rate for l:", accepted.l / i, "\n")

      ### UPDATE sig2 ##########################################################
      sig2new <- rchisq(1, df = sig2old) #rgamma(1, shape = sig2old/2, scale = 2)
      logHR <- logposterior(l = lcurrent, sig2 = sig2new,
                            ystar = ystarold, yprime = yprimeold,
                            given = given) -
        logposterior(l = lcurrent, sig2 = sig2old,
                     ystar = ystarold, yprime = yprimeold,
                     given = given)

      cat("sig2new logHR:", logHR)

      if ( !is.nan(logHR) & !is.na(logHR) ) {

        if (logHR > log(runif(1))) {
          particles.sig2.t[i] <- sig2new
          accepted.sig2 <- accepted.sig2 + 1
        }

      }

      sig2current <- particles.sig2.t[i]
      #cat("\t \t \t -> acceptance rate for sig2:", accepted.sig2 / i, "\n")

      ### UPDATE ystaryprime ###################################################
      # why is this a proposal???
      nugget <- 1e-6  # -7 caused 39/100 warnings
      Rinv <- solve(covMatrix(X = x, sig2 = sig2current,
                              covar.fun = r.matern, l = lcurrent) +
                      diag(nugget, nrow(x)))

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
      # eigen(tau2.xstarprime)$values  # changed nugget to 10e-6 from -8


      ystaryprimenew <- rmvnorm(1, mean = mu.xstarprime, sigma = tau2.xstarprime)
      ystarnew <- ystaryprimenew[1:length(xstar)]
      yprimenew <- ystaryprimenew[(length(xstar)+1):length(ystaryprimenew)]

      particles.ystar.t[i,] <- ystarnew
      particles.yprime.t[i,] <- yprimenew

    }




}



scmc = function(tausteps = seq(from = 0, to = 10, length = 100),
                priorsamples,
                constraint = 0,
                input = suse,
                par = par,
                type = type) {
  #odefn(suse,par,type)
  #timestep=10;prior_fun=odefn;post_fun=m_state_svec;input=suse;par=par;type=type

  theta     = priorsamples
  w         = priorsamples[, 1] * 0 + 1#The value of w   numeric(timestep)
  W         = w / dim(priorsamples)[1] #The value of W
  step = 0
  for (tau in tausteps[-1]) {
    step = step + 1
    w = pnorm(tau * theta) / pnorm(taustep[step + 1] * theta)#updata wj
    W = w * W#update Wj
    # W[,(ts+1)] = (W[,(ts+1)]-mean(W[,(ts+1)]))/sd(W[,(ts+1)]) #Normalize Wj
    ESS = sum(W ^ 2)
    if (ESS < (length(W) / 2)) {
      theta = sample(ori.post, prob = W)#resample
      W     = 1 / dim(priorsamples)[1] #The value of W
    } else{
      theta  = sample(ori.post[, 1], prob = W[, (ts + 1)] / sum(W[, (ts + 1)]))#need to change mh step
    }
  }#end of time step
  after.post = matrix(
    theta[, timestep + 1],
    byrow = F,
    ncol = ncol(post_fun),
    nrow = nrow(post_fun)
  )
  return(after.post)
  b
}
