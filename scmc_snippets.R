## setup
ytrue <- function(x) {log(20*x + 1)}  # x > -1/20; monotone increasing function

given <- list(x = cbind(c(0, 0.1, 0.2, 0.3, 0.4, 0.9, 1)), # inputs where fcn evaluated
              xstar = cbind(seq(0, 1, length.out = 50)), # prediction set inputs
              xprime = cbind(c(0.42, 0.47, 0.52, 0.57, 0.62, # derivative set inputs
                               0.67, 0.72, 0.77, 0.82, 0.87)),
              y = ytrue(c(0, 0.1, 0.2, 0.3, 0.4, 0.9, 1))) # true/known model values

# dealing with list structure
tmp <- list(x = 1:3, y = 3:6, z = 6:9, a = 3:5, b = 3:4, c = 1:10)
tmp
tmp[sample(length(tmp), replace = TRUE)]

getfirst <- unlist(lapply(1:length(tmp), function(x) (tmp[[x]])[1] ))
lapply(getfirst, function(x) x)

# dopar combine argument
library(foreach)
N <- 3
particles.ystar <- lapply(rep(NA, nrow(given$xstar)),
                          function(x) matrix(NA, ncol = 1, nrow = N) )
particles.yprime <- lapply(rep(NA, nrow(given$xprime)),
                           function(x) matrix(NA, ncol = 1, nrow = N) )

combine.particles <- function(out1, out2) {

  return( list(l = c(out1$l, out2$l),
               sig2 = c(out1$sig2, out2$sig2),
               ystar = rbind(out1$ystar, out2$ystar),
               yprime = rbind(out1$yprime, out2$yprime),
               acceptances = rbind(out1$acceptances, out2$acceptances)) )

}

particles.new <- foreach(i = 1:N, .combine = "combine.particles") %dopar% {
  list(l = paste("l:", i), sig2 = paste("sig2:", i),
       ystar = 1:50, yprime = 1:10,
       acceptances = cbind(accepted.l = 1 + sample(0:1),
                           accepted.sig2 = 1 + sample(0:1),
                           accepted.ystaryprime = 1 + sample(0:1)))
}
# str(particles.new)
# particles.new$ystar # each column is the N particles for a single input; ncol = # inputs; nrow = # particles

l.new <- particles.new$l
sig2.new <- particles.new$sig2
ystar.new <- particles.new$ystar
yprime.new <- particles.new$yprime

particles.ystar <- lapply(1:length(given$xstar),
                          function(x) cbind(particles.ystar[[x]], ystar.new[,x]))

particles.yprime <- lapply(1:length(given$xprime),
                           function(x) cbind(particles.yprime[[x]], yprime.new[,x]))


acceptances <- particles.new$acceptances  # N by 3 matrix of acceptances
acceptances
