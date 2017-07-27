# Assessing eta0 results

###
plot(chain.l, type = "l")  # poor mixing for l??
plot(chain.sig2, type = "l")

###
star <- function(i) {
  plot(chain.ystar[i,], type = "l")
  abline(h = chain.ystar[i,1], col = "red")
}

all.equal(chain.ystar[10,], chain.ystar[20,])
star(15)

par(mfrow = c(3, 3))
for (i in 2:10) {star(i)}
par(mfrow = c(1, 1))


###
prime <- function(i) {
  plot(chain.yprime[i,], type = "l")
  abline(h = chain.yprime[i,1], col = "red")
}

prime(1)

par(mfrow = c(3, 3))
for (i in 2:10) {prime(i)}
par(mfrow = c(1, 1))

### Plotting
ymean <- mean(ytrue(x))
curve(log(20*x + 1) - ymean, from = 0, to = 1, lwd = 3)
points(x = x, y = y, pch = 16, col = "green", cex = 3)

for (i in 1:ncol(chain.ystar)) {
  points(x = c(xstar), y = chain.ystar[,i], type = "l", col = "pink")
}

ystar.pred <- apply(chain.ystar, MARGIN = 1, mean)
points(x = seq(0, 1, length.out = 50), y = ystar.pred, type = "l", col = "purple", lwd = 4)
points(x = x, y = y, pch = 16, col = "green", cex = 3)

##
hist(chain.ystar[2,])
abline(v = chain.ystar[2,1], col = "red")


# ### Draws from GP
# ### Set sig2, l values
# l <- mean(chain.l[(length(chain.l) - 10):length(chain.l)])
# sig2 <- mean(chain.sig2)
#
# ymean <- mean(ytrue(x))
# curve(log(20*x + 1) - ymean, from = 0, to = 1, lwd = 3)
# points(x = x, y = y, pch = 16, col = "green", cex = 3)
#
# xstar <- cbind(seq(0, 1, length.out = 60))
# Rinv <- solve(covMatrix(X = x, sig2 = sig2,
#                         covar.fun = r.matern, l = l) +
#                 diag(1e-7, nrow(x)))
#
# R.xstar <- covMatrix(X = xstar, sig2 = sig2, covar.fun = r.matern, l = l) +
#   diag(1e-7, nrow(xstar))
#
#
# # CAREFUL SEE PAPER FOR ARG. ORDER
# r.xstar.x <- covMatrix(X = xstar, X2 = x, sig2 = sig2,
#                                   covar.fun = r.matern, l = l)
# mu.xstar <- r.xstar.x %*% Rinv %*% y
#
# # because covMatrix has sig2 in it, have an extra one multiplying in???
# tau2.xstar <- R.xstar -
#   r.xstar.x %*% Rinv %*% t(r.xstar.x) #/ sig2 #???
#
# # MAKE SURE tau2.xstarprime IS SYMMETRIC
# tau2.xstar <- ( tau2.xstar + t(tau2.xstar) ) / 2
#
# # MAKE SURE tau2.xstarprime IS POSTIVE SEMI-DEFINITE
# # eigen(tau2.xstarprime)$values  # changed nugget to 10e-6 from -8
#
# posterior.draws <- rmvnorm(10, mean = mu.xstar, sigma = tau2.xstar)
# all.equal(posterior.draws[1,], posterior.draws[2,])
#
# for (i in 1:nrow(posterior.draws)) {
#   points(x = c(xstar), y = posterior.draws[i,], type = "l", col = "pink", lwd = 2)
# }
#
# points(x = x, y = y, pch = 16, col = "green", cex = 3)
