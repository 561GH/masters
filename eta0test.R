## run setup

# monotone increasing function
ytrue <- function(x) {log(20*x + 1)}  # x > -1/20

# inputs where fcn evaluated
x <- c(0, 0.1, 0.2, 0.3, 0.4, 0.9, 1)
# true/known model values
y <- ytrue(x)

# prediction set inputs
xstar <- seq(0, 1, length.out = 50)

# derivative set inputs
xprime <- c(0.42, 0.47, 0.52, 0.57, 0.62, 0.67, 0.72, 0.77, 0.82, 0.87)

given <- list(x = cbind(x), xprime = cbind(xprime), xstar = cbind(xstar),
              y = y)

## start MH within gibbs

eta.init <- list(l = 5.7, sig2 = 1500,
                 ystar = ytrue(xstar),
                 yprime = 20 / (20 * xprime))

N <- 100
v1 <- 1
v2 <- 1
chain.l <- rep(NA, N)
chain.sig2 <- rep(NA, N)
chain.ystar <- matrix(NA, nrow = length(xstar), ncol = N)
chain.yprime <- matrix(NA, nrow = length(xprime), ncol = N)

for (i in 1:N) {

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
  logHR <- log.posterior(l = lnew,
                      sig2 = sig2old, ystar = ystarold, yprime = yprimeold,
                      given = given) -
    log.posterior(l = lold,
                  sig2 = sig2old, ystar = ystarold, yprime = yprimeold,
                  given = given)

  cat("lnew logHR:", logHR, "\n")

  if (logHR > log(runif(1))) {
    chain.l[i] <- lnew
  } else {
    chain.l[i] <- lold
  }

  lcurrent <- chain.l[i]

  ## UPDATE sig2 ###############################################################
  sig2new <- rgamma(1, shape = sig2old/2, scale = 2)
  logHR <- log.posterior(l = lcurrent, sig2 = sig2new,
                         ystar = ystarold, yprime = yprimeold,
                         given = given) -
    log.posterior(l = lcurrent, sig2 = sig2old,
                  ystar = ystarold, yprime = yprimeold,
                  given = given)

  cat("sig2new logHR:", logHR, "\n")

  if (logHR > log(runif(1))) {
    chain.sig2[i] <- sig2new
  } else {
    chain.sig2[i] <- sig2old
  }

  sig2current <- chain.sig2[i]

  ## UPDATE ystaryprime ########################################################


  ystaryprimenew <- rmvnorm(1, mean = mystatyprime, sigma = v2 * S)
  logHR <- log.posterior(l = lcurrent, sig2 = sig2new,
                         ystar = ystarold, yprime = yprimeold,
                         given = given) -
    log.posterior(l = lcurrent, sig2 = sig2old,
                  ystar = ystarold, yprime = yprimeold,
                  given = given)

  cat("sig2new logHR:", logHR, "\n")

  if (logHR > log(runif(1))) {
    chain.sig2[i] <- sig2new
  } else {
    chain.sig2[i] <- sig2old
  }

  sig2current <- chain.sig2[i]

}
