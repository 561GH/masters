# In parallel
clusters <- parallel::detectCores() / 2
cl <- parallel::makeCluster(clusters)
doSNOW::registerDoSNOW(cl)
## run stuff here
parallel::stopCluster(cl)


# Plotting results from scmc_ex1
ystar.true <- ytrue(given$xstar) - mean(ytrue(given$x))

plot(x = given$x, y = given$y, pch = 16, cex = 2)
for (n in 1:N) {
  # points(x = given$xstar, y = particles.ystar[n,14,], type = "l", col = "red")
  # points(x = given$xstar, y = particles.ystar[n,13,], type = "l", col = "green")
  # points(x = given$xstar, y = particles.ystar[n,12,], type = "l", col = "blue")
  points(x = given$xstar, y = particles.ystar[n,14,], type = "l", col = "red")
}
points(x = given$xstar, y = ytrue(given$xstar) - mean(ytrue(given$x)), type = "l")




plot(x = given$xprime, y = 20 / (20 * given$xprime), pch = 16, cex = 2)
for (n in 1:N) {
  points(x = given$xprime, y = particles.yprime[n,7,], type = "l", col = n)
}

tmp <- apply(particles.yprime[,13,], MARGIN = 2, FUN = mean)
plot(x = given$xprime, y = tmp, type = "l")
tmp
points(x = given$xprime, y = c(20 / (20 * given$xprime)))

################################################################################
library(MonGP)
library(MASS)
library(lgarch)
test <- SMC_MonGP(x = given$x, y_obs = given$y, xd = given$xprime, newx = given$xstar,
                  initial_values = list(l = 5.7, sig2 = 1), N = 5, L = 20)
