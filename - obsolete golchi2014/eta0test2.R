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

given <- list(x = cbind(x), xprime = cbind(xprime), xstar = cbind(xstar),
              y = y)

## start MH within gibbs #######################################################

eta.init <- list(l = 0.3, sig2 = 1, #l = 5.76, sig2 = 5.76,
                 ystar = ytrue(xstar) - mean(ytrue(x)),
                 yprime = 20 / (20 * xprime) )

set.seed(561)
N <- 100  # iterations
v1 <- 8  # step size for l proposal
v2 <- 8  # step size for sig2 proposal
chain.l <- rep(NA, N)
chain.sig2 <- rep(NA, N)
chain.ystar <- matrix(NA, nrow = length(xstar), ncol = N)
chain.yprime <- matrix(NA, nrow = length(xprime), ncol = N)

### Track acceptance rate
accepted.l <- accepted.sig2 <- 0

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
  cat("    -> acceptance rate for l:", accepted.l / i, "\n")

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
  cat("    -> acceptance rate for sig2:", accepted.sig2 / i, "\n")

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
#)
# saveRDS(out, file = "/Users/ggh/Git/masters/eta0out.rds")

# after MH within Gibbs
chains <- list(eta.init = eta.init,
               chain.l = chain.l,
               chain.sig2 = chain.sig2,
               chain.ystar = chain.ystar,
               chain.yprime = chain.yprime)
#saveRDS(chains, file = "/Users/ggh/Git/masters/eta0chains.rds")

#
# chain.l
# chain.sig2
# cbind(eta.init$yprime, chain.yprime[,ncol(chain.yprime)])
#
# head(cbind(eta.init$ystar, chain.ystar))
# pred.ystar <- apply(chain.ystar[,(ncol(chain.ystar) - 100):ncol(chain.ystar)],
#                     MARGIN = 1, FUN = mean)
# cbind(ytrue(xstar) - mean(ytrue(x)), pred.ystar - mean(pred.ystar))
#
#
# # draw
# plot(x = x, y = ytrue(x) - mean(ytrue(x)), pch = 16)
# points(x = xstar, y = ytrue(xstar)-mean(ytrue(x)), pch = 16, col = "green")
#
# R <- covMatrix(X = x, sig2 = mean(chain.sig2),
#                covar.fun = r.matern, l = mean(chain.l)) + diag(rep(0.0000001, 7))
# r.xstar.x <- covMatrix(X = xstar, X2 = x, sig2 = mean(chain.sig2),
#                        covar.fun = r.matern, l = mean(chain.l))
# mu.star <- r.xstar.x  %*% solve(R) %*% y
# cov.star <- mean(chain.sig2) * ( 1 -  r.xstar.x %*% solve(R) %*% t(r.xstar.x) ) #+ diag(rep(40, 50))
#
# eig <- eigen(cov.star)
# evs <- eig$values
# evs
# all.equal(cov.star, t(cov.star))
# pred.draws <- rmvnorm(5, mean = mu.star, sigma = cov.star)
#
# for (i in 1:nrow(pred.draws)) {
#   points(x = xstar, y = pred.draws[i,], col = "red", type = "l")
# }
# ###
#
# Kystary <- covMatrix(X = xstar, X2 = x,
#                      covar.fun = r.matern, sig2 = 1, l = 1)
# Kyy <- covMatrix(X = x,
#                  covar.fun = r.matern, sig2 = 1, l = 1)
# Kystarystar <- covMatrix(X = xstar,
#                          covar.fun = r.matern, sig2 = 1, l = 1)
#
# pred.mean <- Kystary %*% solve(Kyy) %*% (y - 0)
# points(x = xstar, y = pred.mean, col = "red", cex = 2)
#
# small.term <- Kystary %*% solve(Kyy) %*% t(Kystary)
# pred.covar <- Kystarystar - small.term + 0.01
# all.equal(pred.covar, t(pred.covar))
# pred.draws <- rmvnorm(5, mean = pred.mean, sigma = pred.covar)
#
# for (i in 1:nrow(pred.draws)) {
#   points(x = xstar, y = pred.draws[i,], col = "red", type = "l")
# }
#
#
# S11 <- covMatrix(X = xstar, sig2 = mean(chain.sig2),
#                  covar.fun = r.matern, l = mean(chain.l))
# S22 <- covMatrix(X = xprime, sig2 = mean(chain.sig2),
#                  covar.fun = r.matern2, l = mean(chain.l))
# # CAREFUL SEE PAPER FOR ARG. ORDER
# S21 <- covMatrix(X = xprime, X2 = xstar, sig2 = mean(chain.sig2),
#                  covar.fun = r.matern1, l = mean(chain.l))
# S <- rbind(cbind(S11, t(S21)),
#            cbind(S21, S22))
#
# draws.ystar <- rmvnorm(5, mean = pred.ystar - mean(pred.ystar), sigma = S11)
#
# plot(x = x, y = ytrue(x) - mean(ytrue(x)))
# points(x = xstar, y = ytrue(xstar) - mean(ytrue(x)), col = "green")
# #points(x = xstar, y = pred.ystar - mean(pred.ystar))
#
# for (i in 1:nrow(draws.ystar)) {
#   points(x = xstar, y = draws.ystar[i,], pch = 16)
# }
# points(x = xstar, y = pred.ystar - mean(pred.ystar), pch = 16, col = "red")
#
#
