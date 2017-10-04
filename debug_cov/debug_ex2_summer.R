# debugging ex2

## testing columns of tau2.xstar
cos.angle <- matrix(rep(NA, ncol(tau2.xstarprime)^2), ncol = ncol(tau2.xstarprime))
for (i in 1:(ncol(tau2.xstarprime) - 1)) {
  for (j in (i+1):ncol(tau2.xstarprime)) {
    cos.angle[i, j] <- sum(tau2.xstarprime[,i] * tau2.xstarprime[,j]) /
      ( sqrt( sum(tau2.xstarprime[,i]^2) ) * sqrt( sum(tau2.xstarprime[,j]^2) ) )
  }
}
hist(cos.angle)
sort(cos.angle)
which(cos.angle > 0.98, arr.ind = TRUE)

# ## testing with MCMC_unconstrained
# x <- given$x
# y_obs <- given$y
# xd <- given$xprime
# newx <- given$xstar
# N <- 10
# burn <- 0
# initial_values <- initial.values

## following Wang lemma 3.1 ####################################################

# given: x, xstar, xprime, y; y must be drawn from 0 mean GP
# l: length-scale GP parameter
# sig2: constant variance GP parameter

x <- given$x
xstar <- given$xstar
xprime <- given$xprime
y <- given$y

R <- covMatrix(X1 = x, X2 = x, sig2 = sig2, l = l,
               covar.fun = "matern") +
  diag(nugget1, nrow(x))
Linv <- solve( t(chol(R)) )  # recall need transpose to match std chol def
Rinv <- t(Linv) %*% Linv

A2 <- covMatrix(X1 = xprime1, X2 = xstar,
                sig2 = sig2, l = l, covar.fun = "matern1", d1 = 1)

m1 <-


n1 <- nrow(xprime) / 2  # number of derivs. in d1
xprime1 <- xprime[1:n1,]
xprime2 <- xprime[-(1:n1),]

S11 <- covMatrix(X1 = xstar, X2 = xstar,
                 sig2 = sig2, l = l, covar.fun = "matern")
S22 <- covMatrix(X1 = xprime1, X2 = xprime1,
                 sig2 = sig2, l = l, covar.fun = "matern2", d1 = 1, d2 = 1)
S33 <- covMatrix(X1 = xprime2, X2 = xprime2,
                 sig2 = sig2, l = l, covar.fun = "matern2", d1 = 2, d2 = 2)

S21 <- covMatrix(X1 = xprime1, X2 = xstar, # CAREFUL SEE PAPER FOR ARG. ORDER
                 sig2 = sig2, l = l, covar.fun = "matern1", d1 = 1)
S31 <- covMatrix(X1 = xprime2, X2 = xstar, # CAREFUL SEE PAPER FOR ARG. ORDER
                 sig2 = sig2, l = l, covar.fun = "matern1", d1 = 2)

S32 <- covMatrix(X1 = xprime2, X2 = xprime1,
                 sig2 = sig2, l = l, covar.fun = "matern2", d1 = 2, d2 = 1)

R.xstarxprime <- rbind(cbind(S11, t(S21), t(S31)),
                       cbind(S21,   S22,  t(S32)),
                       cbind(S31,   S32,    S33)) +
  diag(rep(nugget2, nrow(xstar) + nrow(xprime)))

# CAREFUL SEE PAPER FOR ARG. ORDER
r.xstarprime.x <- rbind(covMatrix(X1 = xstar, X2 = x,
                                  sig2 = sig2, l = l,
                                  covar.fun = "matern"),
                        covMatrix(X1 = xprime1, X2 = x,
                                  sig2 = sig2, l = l,
                                  covar.fun = "matern1", d1 = 1),
                        covMatrix(X1 = xprime2, X2 = x,
                                  sig2 = sig2, l = l,
                                  covar.fun = "matern1", d1 = 2))
mu.xstarprime <- r.xstarprime.x %*% Rinv %*% y

# Following Rasmussen (2.19)
tau2.xstarprime <- R.xstarxprime -
  r.xstarprime.x %*% Rinv %*% t(r.xstarprime.x)

# MAKE SURE tau2.xstarprime IS SYMMETRIC
tau2.xstarprime <- ( tau2.xstarprime + t(tau2.xstarprime) ) / 2


return(list(mu.xstarprime = mu.xstarprime,
            cov.xstarprime = tau2.xstarprime,
            R.xx = R))

