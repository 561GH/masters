## run setup
library(mvtnorm)

# monotone increasing function
ytrue <- function(x) {log(20*x + 1)}  # x > -1/20

# inputs where fcn evaluated
x <- cbind(c(0, 0.1, 0.2, 0.3, 0.4, 0.9, 1))
# true/known model values with 0 mean
y <- ytrue(x) - mean(ytrue(x))

# prediction set inputs
xstar <- cbind(seq(0, 1, length.out = 50))

# derivative set inputs
xprime <- cbind(c(0.42, 0.47, 0.52, 0.57, 0.62, 0.67, 0.72, 0.77, 0.82, 0.87))

## start MH within gibbs #######################################################

eta0 <- function(eta.init = init,
                 given = list(x = cbind(x),
                              xprime = cbind(xprime), xstar = cbind(xstar),
                              y = y),  # data, locations, obs, etc.
                 N = 100, # iterations
                 v1 = 8,  # step size for l proposal
                 v2 = 8) {# step size for sig2 proposal

  chain.l <- rep(NA, N)
  chain.sig2 <- rep(NA, N)
  chain.ystar <- matrix(NA, nrow = length(xstar), ncol = N)
  chain.yprime <- matrix(NA, nrow = length(xprime), ncol = N)

  ### Track acceptance rate
  accepted.l <- accepted.sig2 <- 0

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

    chain.l[i] <- lold
    chain.sig2[i] <- sig2old
    chain.ystar[,i] <- ystarold
    chain.yprime[,i] <- yprimeold

    ## UPDATE l ##################################################################
    lnew <- rnorm(1, mean = lold, sd = v1)

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

    cat("sig2new logHR:", logHR)

    if ( !is.nan(logHR) & !is.na(logHR) ) {

      if (logHR > log(runif(1))) {
        chain.sig2[i] <- sig2new
        accepted.sig2 <- accepted.sig2 + 1
      }

    }

    sig2current <- chain.sig2[i]
    cat("\t \t \t -> acceptance rate for sig2:", accepted.sig2 / i, "\n")

    ## UPDATE ystaryprime ########################################################
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

    chain.ystar[,i] <- ystarnew
    chain.yprime[,i] <- yprimenew

    cat("Finished iteration:", i, "\n \n")

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

## tuning ######################################################################
init.lsig2 <- c(5.7, 1)
init <- list(l = 5.7, sig2 = 1,
             ystar = ytrue(xstar) - mean(ytrue(x)),
             yprime = 20 / (20 * xprime) )

trials.v1 <- 0.01 * 2^(-2:2)
trials.v2 <- 6 * 2^(-2:2)
trials.v1v2 <- expand.grid(trials.v1, trials.v2)


steps <- 50
runs <- NULL
for (j in 1:14) {

  rate.l <- rate.sig2 <- NULL
  for (i in 1:nrow(trials.v1v2)) {
    test <- eta0(v1 = trials.v1v2[i, 1], v2 = trials.v1v2[i, 2], N = steps)
    rate.l <- c(rate.l, test$accepted.l / steps)
    rate.sig2 <- c(rate.sig2, test$accepted.sig2 / steps)
    cat("Finished test:", i, "of round:", j, "\n")
  }

  runs <- rbind(runs, cbind(1:25, rate.l, rate.sig2))
}

runs <- data.frame(runs)
colnames(runs) <- c("combo.row", "rate.l", "rate.sig2")

res10 <- rbind(cbind(combo.row = 1:25, rate.l = res$rate.l, rate.sig2 = res$rate.sig2), runs)
res10$combo.row <- factor(res10$combo.row)
#write.csv(res10, file = "/Users/ggh/Desktop/res10.csv", row.names = FALSE)

# plots
library(ggplot2)

## 10 rounds
res10gg <- tidyr::gather(data = res10, parameter, rate, -1)
gg <- ggplot(data = res10gg, aes(x = combo.row, y = rate))
gg + geom_boxplot(aes(fill = parameter)) + geom_hline(yintercept = 0.4)

gg + geom_boxplot() + facet_wrap(~ parameter) + geom_hline(yintercept = 0.4)

trials.v1v2[8,]
#   Var1   Var2
# 8 0.01    3
#' - somehwere in between here for l step size
#' - for sig2, steps tested seem to affect spread more than anything?

# Generate the particles for SCMC
