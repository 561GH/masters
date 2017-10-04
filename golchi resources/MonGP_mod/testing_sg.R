library(MASS); library(lgarch); library(dplyr); library(ggplot2)
ytrue <- function(x) {log(20 * x + 1)}  # x > -1/20; monotone increasing function
given <- list(x = cbind(c(0, 0.1, 0.2, 0.3, 0.4, 0.9, 1)),
              xstar = cbind(seq(0, 1, length.out = 50)),
              xprime = cbind(c(0.42, 0.47, 0.52, 0.57, 0.62,
                               0.67, 0.72, 0.77, 0.82, 0.87)),
              y = ytrue(c(0, 0.1, 0.2, 0.3, 0.4, 0.9, 1)) -
                mean(ytrue(c(0, 0.1, 0.2, 0.3, 0.4, 0.9, 1))) )

x <- given$x
newx <- given$xstar
y_obs <- given$y
xd <- given$xprime

obj <- SMC_MonGP(x, y_obs, xd, newx,
                 N = 100, L = 10, list(l = .5, sig2 = 10))
sample_path(obj)
debug <- obj$debug
W <- obj$weights
boxplot(W)
particles.yprime <- obj$yd
for (t in 2:10) {
  print( sum(particles.yprime[,t,] %in% particles.yprime[,1,]) )
}

particles.ystar <- obj$y
for (t in 2:10) {
  print( sum(particles.ystar[,t,] %in% particles.ystar[,10,]) )
}

# seems like for each time step, the only important part is the resampling????

  l <- obj$l
  sig2 <- obj$sig2

  hist(l[,10])
  hist(sig2[,10])

  med.l <- median(l[,10])
  med.sig2 <- median(sig2[,10])
