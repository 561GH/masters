pred.paras <- predictiveMeanCov(given = given, l = 0.5, sig2 = 1)
cov <- pred.paras$cov.xstarprime
diag(cov)
options(digits = 2)

x <- given$x
xstar <- given$xstar
xprime <- given$xprime
y <- given$y
sig2 <- 1

nugget <- 1e-6
R <- covMatrix(X1 = x, X2 = x, sig2 = sig2, l = l,
               covar.fun = "matern")
Linv <- solve( t(chol(R)) )  # recall need transpose to match std chol def
Rinv <- t(Linv) %*% Linv

S11 <- covMatrix(X1 = xstar, X2 = xstar, sig2 = sig2, l = l,
                 covar.fun = "matern")
S22 <- covMatrix(X1 = xprime, X2 = xprime, sig2 = sig2, l = l,
                 covar.fun = "matern2", d1 = 1, d2 = 1)
S21 <- covMatrix(X1 = xprime, X2 = xstar,  # CAREFUL SEE PAPER FOR ARG. ORDER
                 sig2 = sig2, l = l,
                 covar.fun = "matern1", d1 = 1)
R.xstarxprime <- rbind(cbind(S11, t(S21)),
                       cbind(S21, S22))

r.xstarprime.x <- rbind(covMatrix(X1 = xstar, X2 = x,
                                  sig2 = sig2, l = l,
                                  covar.fun = "matern"),
                        covMatrix(X1 = xprime, X2 = x,
                                  sig2 = sig2, l = l,
                                  covar.fun = "matern1", d1 = 1))

# because covMatrix has sig2 in it, have an extra one multiplying in???
tau2.xstarprime <- R.xstarxprime -
  r.xstarprime.x %*% Rinv %*% t(r.xstarprime.x) #/ sig2current #???

V <- diag(diag(tau2.xstarprime))
Vsqrt.inv <- solve( (sqrt(V)) )
corr.tau2.xstarprime <- Vsqrt.inv %*% tau2.xstarprime %*% Vsqrt.inv

###
x <- given$x
newx <- given$xstar
xd <- given$xprime
l <- 0.5
library(MASS)

R = matcor(x, x, l)
Rinv = ginv(R)
c0 = matcor(x, newx, l)
c1 = matcorder(x, xd, l, 1)
C = cbind(c0, c1)
c00 = matcor(newx, newx, l)
c01 = matcorder(newx, xd, l, 1)
c11 = matsecder(xd, xd, l, 1, 1)
M = rbind(cbind(c00, c01), cbind(t(c01), c11))
corrmat = M - t(C) %*% Rinv %*% C

all.equal(corrmat, tau2.xstarprime)
