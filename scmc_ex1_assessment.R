# In parallel
# clusters <- parallel::detectCores() / 2
# cl <- parallel::makeCluster(clusters)
# doSNOW::registerDoSNOW(cl)
# ## run stuff here
# parallel::stopCluster(cl)

library(ggplot2); library(dplyr)
plotSamplePaths <- function(given,
                            particles.ystar) {
  M <- dim(particles.ystar)[2]
  N <- dim(particles.ystar)[1]
  dimnames(particles.ystar) <- list(sample = 1:N, step = 1:M, point = given$xstar)
  df0 <- as.data.frame.table(particles.ystar)
  names(df0)[4] <- "y"
  select_step <- floor(seq(1, M, length = 5))
  df0 <- df0 %>% filter(step %in% select_step) %>%
    group_by(point, step) %>%
    summarise(mean = mean(y),
              lower = quantile(y, .025),
              upper = quantile(y, 0.975))
  df0 <- as.data.frame(df0)

  ystar.true <- ytrue(given$xstar) - mean(ytrue(given$x))
  truth <- data.frame(point = factor(given$xstar),
                      step = factor(rep("truth", nrow(given$xstar))),
                      mean = ystar.true,
                      lower = ystar.true,
                      upper = ystar.true)

  observed <- data.frame(point = factor(given$x),
                         step = factor(rep("obs", nrow(given$x))),
                         mean = given$y,
                         lower = given$y,
                         upper = given$y)

  # bind truth to df0 so that colours steps in the desired order
  #df0 <- rbind(df0, truth)
  df0$step <- droplevels(df0$step)

  p <- ggplot(df0, aes(x = as.numeric(as.character(point)),
                      y = mean, fill = step)) +
    geom_line() +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3) +
    geom_line(colour = "darkRed", data = truth, size = 0.6) +
    geom_point(aes(x = as.numeric(as.character(point)), y = mean),
               data = observed, color = "darkRed") +
    scale_fill_brewer(palette = "YlOrBr") +
    xlab("x") + ylab("y")

  return(p)
}


# Plotting results from scmc_ex1
t <- M
plotSamplePaths(given = given,
                particles.ystar = particles.ystar[,1:t,])

# step.sizes
# tail(acceptances[[2]])
# tail(acceptances[[floor(M/2)]])
# tail(acceptances[[M]])
#
# hist(particles.l[,t])
# hist(particles.sig2[,t])
#
# median(particles.l[,t])
# median(particles.sig2[,t])

#ystar.true <- ytrue(given$xstar) - mean(ytrue(given$x))

plot(x = given$x, y = given$y, pch = 16, cex = 2) #, ylim = c(-0.4, 0.7))
for (n in 1:N) {
  points(x = given$xstar, y = particles.ystar[n,M,], type = "l", col = "red")
  points(x = given$xstar, y = particles.ystar[n,10,], type = "l", col = "green")
  points(x = given$xstar, y = particles.ystar[n,5,], type = "l", col = "blue")
  points(x = given$xstar, y = particles.ystar[n,1,], type = "l", col = "pink")
}
#points(x = given$xstar, y = ytrue(given$xstar) - mean(ytrue(given$x)), type = "l")


# plot(x = given$xprime, y = yprimetrue(given$xprime), pch = 16, cex = 2)
# for (n in 1:N) {
#   points(x = given$xprime, y = particles.yprime[n,8,], col = "blue")
# }
# points(x = given$xprime, y = apply(particles.yprime[,5,], 2, mean), col = "red", cex = 2, pch = 16)
# points(x = given$xprime, y = apply(particles.yprime[,10,], 2, mean), col = "green", cex = 2, pch = 16)
# points(x = given$xprime, y = apply(particles.yprime[,14,], 2, mean), col = "blue", cex = 2, pch = 16)
# points(x = given$xprime, y = apply(particles.yprime[,t,], 2, mean), col = "pink", cex = 2, pch = 16)
#
# points(x = xd, y = apply(yd[,5,], MARGIN = 2, FUN = mean), col = "red", cex = 2, pch = 17)
# points(x = xd, y = apply(yd[,10,], MARGIN = 2, FUN = mean), col = "green", cex = 2, pch = 17)
#
# points(x = given$xprime, y = c(20 / (20 * given$xprime)))

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

library(MonGP)
obj <- readRDS("/Users/ggh/Desktop/sg.rds")
sample_path(obj)

l <- obj$l
hist(l[,10])
sig2 <- obj$sig2
hist(sig2[,1])

# more tests on shirin's code
library(MonGP); library(MASS); library(lgarch); library(dplyr); library(ggplot2)
x <- given$x
newx <- given$xstar
y_obs <- given$y
xd <- given$xprime

med.l <- med.sig2 <- rep(NA, 10)
for (i in 1:10) {
  obj <- SMC_MonGP(x, y_obs, xd, newx,
                   N = 100, L = 10, list(l = .5, sig2 = 10))
  l <- obj$l
  sig2 <- obj$sig2

  pdf(paste("/Users/ggh/Desktop/sg/hist", i, ".pdf", sep = ""), width = 14)
  par(mfrow = c(1, 2))
  hist(l[,10])
  hist(sig2[,10])
  dev.off()

  ggsave(plot = sample_path(obj),
         filename = paste("/Users/ggh/Desktop/sg/plot", i, ".pdf", sep = ""))

  med.l[i] <- median(l[,10])
  med.sig2[i] <- median(sig2[,10])
}
saveRDS(cbind(med.l, med.sig2), file = "/Users/ggh/Desktop/sg/medians.rds")
cbind(med.l, med.sig2)
# so the code above will randomly die sometimes; not sure why

# compare with mine T_T
i <- 1
med.l.gh <- med.sig2.gh <- rep(NA, 10)

###
i <- i + 1
pdf(paste("/Users/ggh/Desktop/sg/hist_gh", i, ".pdf", sep = ""), width = 14)
par(mfrow = c(1, 2))
hist(particles.l[,t])
hist(particles.sig2[,t])
dev.off()


ggsave(plot = plotSamplePaths(given = given,
                              particles.l = particles.l[,1:t],
                              particles.ystar = particles.ystar[,1:t,]),
       filename = paste("/Users/ggh/Desktop/sg/plot_gh", i, ".pdf", sep = ""))

med.l.gh[i] <- median(particles.l[,t])
med.sig2.gh[i] <- median(particles.sig2[,t])


saveRDS(cbind(med.l.gh, med.sig2.gh), file = "/Users/ggh/Desktop/sg/medians_gh.rds")
params.gh <- cbind(med.l.gh, med.sig2.gh)
params.sg <- cbind(med.l, med.sig2)

boxplot(params.gh)
boxplot(params.sg)

boxl <- cbind(med.l, med.l.gh)
boxsig2 <- cbind(med.sig2, med.sig2.gh)

boxplot(boxl)
boxplot(boxsig2)

plot(med.sig2, ylim = c(0, 19))
points(med.l, col = "blue")
points(med.sig2.gh, pch = 16)
points(med.l.gh, pch = 16, col = "blue")
