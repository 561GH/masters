# In parallel
clusters <- parallel::detectCores() / 2
cl <- parallel::makeCluster(clusters)
doSNOW::registerDoSNOW(cl)
## run stuff here
parallel::stopCluster(cl)


# Plotting results from scmc_ex1
library(ggplot2); library(dplyr)
plotSamplePaths(given = given,
                particles.l = particles.l[,1:t],
                particles.ystar = particles.ystar[,1:t,])

hist(particles.l[,t])
hist(particles.sig2[,t])

ystar.true <- ytrue(given$xstar) - mean(ytrue(given$x))

# plot(x = given$x, y = given$y, pch = 16, cex = 2)
# for (n in 1:N) {
#   #points(x = given$xstar, y = particles.ystar[n,10,], type = "l", col = "red")
#   points(x = given$xstar, y = particles.ystar[n,13,], type = "l", col = "green")
#   points(x = given$xstar, y = particles.ystar[n,14,], type = "l", col = "blue")
#   points(x = given$xstar, y = particles.ystar[n,19,], type = "l", col = "pink")
# }
# points(x = given$xstar, y = ytrue(given$xstar) - mean(ytrue(given$x)), type = "l")


plot(x = given$xprime, y = 20 / (20 * given$xprime), pch = 16, cex = 2)
# for (n in 1:N) {
#   points(x = given$xprime, y = particles.yprime[n,14,], col = n)
# }
points(x = given$xprime, y = apply(particles.yprime[,5,], 2, mean), col = "red", cex = 2, pch = 16)
points(x = given$xprime, y = apply(particles.yprime[,10,], 2, mean), col = "green", cex = 2, pch = 16)
points(x = given$xprime, y = apply(particles.yprime[,14,], 2, mean), col = "blue", cex = 2, pch = 16)
points(x = given$xprime, y = apply(particles.yprime[,t,], 2, mean), col = "pink", cex = 2, pch = 16)

points(x = xd, y = apply(yd[,5,], MARGIN = 2, FUN = mean), col = "red", cex = 2, pch = 17)
points(x = xd, y = apply(yd[,10,], MARGIN = 2, FUN = mean), col = "green", cex = 2, pch = 17)

points(x = given$xprime, y = c(20 / (20 * given$xprime)))

################################################################################
# library(MonGP); library(MASS); library(lgarch); library(dplyr); library(ggplot2)
# x <- as.matrix(c(seq(0,.4,.1), seq(.9,1,.1)))
# fnct <- function(x) {
#   out <- log(20 * x + 1)
#   return(out)
# }
# data <- fnct(x)
# mu <- mean(data)
# se <- as.numeric(sqrt(var(data)))
# y_obs <- (data - mu)/se
# newx <- as.matrix(seq(.01,.99,length=50))
# xd <- as.matrix(seq(.42,.87,length=10))
# initial_values = list(l = .5, sig2 = 10)
# obj = SMC_MonGP(x, y_obs, xd, newx, N = 1000, L = 10, initial_values)
# sample_path(obj)

obj <- readRDS("/Users/ggh/Desktop/sg.rds")
l <- obj$l
hist(l[,10])
sig2 <- obj$sig2
hist(sig2[,1])

xd <- as.matrix(seq(.42,.87,length=10))
yd <- obj$yd
ystar <- obj$y

plot(x = xd, y = 20 / (20 * xd))
points(x = xd, y = apply(yd[,10,], MARGIN = 2, FUN = mean), col = "red")

sample_path(obj)


