# Setup EX2

ytrue <- function(x1, x2) {
  11*x1^10 + 9*x1^8 + 7*x1^6 +
    10*x2^9 + 8*x2^7
}

yprime1 <- function(x1, x2) {11*10*x1^9 + 9*8*x1^7 + 7*6*x1^5}
yprime2 <- function(x1, x2) {10*9*x2^8 + 8*7*x2^6}

# Plots in 3D ##################################################################
# library(plot3D)
# curve_3d <- function(f2,
#                      x_range = c(0, 1),
#                      y_range = c(0, 1),
#                      col = 1:6) {
#
#   x1 <- seq(x_range[1], x_range[2], len = 30)
#   x2 <- seq(y_range[1], y_range[2], len = 30)
#   y <- outer(x1, x2, FUN = f2)
#
#   # persp3D(xvec, yvec, fz,
#   #         contour = TRUE, theta = 15, phi = 40,
#   #         axes = TRUE, scale = 2, box = TRUE, nticks = 5,
#   #         ticktype = "detailed")
#   #
#   persp(x1, x2, y, theta = 15, phi = 40,
#         shade = 0.4, col = "cyan",
#         axes = TRUE, scale = 2, box = TRUE, nticks = 5,
#         ticktype = "detailed")
#
# }
#
# curve_3d(ytrue)
# curve_3d(yprime1)
# curve_3d(yprime2)

# neex xquartz?
# require(rgl)
#
# points3d(0.5,0.5,20,col="red")

# # Setup ########################################################################
# library(DiceDesign)
#
# x <- lhsDesign(n = 15, dimension = 2, seed = 761)$design
# y <- ytrue(x1 = x[,1], x2 = x[,2])
# plot(x = x[,1], y = x[,2], pch = 16)
#
# xstar <- rbind(c(3, 5), c(4, 3), c(5, 4), c(6, 7), c(7, 6)) / 10
# points(x = xstar[,1], y = xstar[,2], pch = 16, col = "blue")
#
# get.xprime1 <- function(xstar) {
#   rbind(cbind(xstar[1] - 0.03, xstar[2]),
#         cbind(xstar[1] - 0.06, xstar[2]),
#         cbind(xstar[1] + 0.03, xstar[2]),
#         cbind(xstar[1] + 0.06, xstar[2]))
# }
# get.xprime2 <- function(xstar) {
#   rbind(cbind(xstar[1], xstar[2] - 0.03),
#         cbind(xstar[1], xstar[2] - 0.06),
#         cbind(xstar[1], xstar[2] + 0.03),
#         cbind(xstar[1], xstar[2] + 0.06))
# }
#
# xprime1 <- xprime2 <- yprime.init1 <- yprime.init2 <- NULL
# for (i in 1:nrow(xstar)) {
#   tmp1 <- get.xprime1(xstar[i,])
#   tmp2 <- get.xprime2(xstar[i,])
#   xprime1 <- rbind(xprime1, tmp1)
#   xprime2 <- rbind(xprime2, tmp2)
#   yprime.init1 <- c(yprime.init1, yprime1(tmp1[, 1], tmp1[, 2]))
#   yprime.init2 <- c(yprime.init2, yprime2(tmp2[, 1], tmp2[, 2]))
# }
# points(x = xprime1[,1], y = xprime1[,2], col = "red", pch = 16)
# points(x = xprime2[,1], y = xprime2[,2], col = "red", pch = 16)
#
# given <- list(x = x,
#               xstar = xstar,
#               xprime = rbind(xprime1, xprime2),
#               y = y - mean(y))
#
# initial.values <- list(l = c(0.5, 0.5), sig2 = 10,
#                        ystar = ytrue(xstar[,1], xstar[,2]) - mean(y) +
#                          rnorm(nrow(xstar), 0, 0.5),
#                        yprime = c(yprime.init1, yprime.init2))
#
# rm(x, y, tmp1, tmp2, xprime1, xprime2, xstar,
#    yprime.init1, yprime.init2, i, get.xprime1, get.xprime2)

# Setup reduced ################################################################
library(DiceDesign)

x <- lhsDesign(n = 15, dimension = 2, seed = 761)$design
y <- ytrue(x1 = x[,1], x2 = x[,2])
plot(x = x[,1], y = x[,2], pch = 16)

xstar <- rbind(c(3, 5), c(4, 3), c(5, 4), c(6, 7), c(7, 6)) / 10
points(x = xstar[,1], y = xstar[,2], pch = 16, col = "blue")

get.xprime1 <- function(xstar) {
  rbind(cbind(xstar[1] - 0.03, xstar[2]),
        cbind(xstar[1] - 0.06, xstar[2]),
        cbind(xstar[1] + 0.03, xstar[2]),
        cbind(xstar[1] + 0.06, xstar[2]))
}
get.xprime2 <- function(xstar) {
  rbind(cbind(xstar[1], xstar[2] - 0.03),
        cbind(xstar[1], xstar[2] - 0.06),
        cbind(xstar[1], xstar[2] + 0.03),
        cbind(xstar[1], xstar[2] + 0.06))
}

xprime1 <- xprime2 <- yprime.init1 <- yprime.init2 <- NULL
for (i in 1:nrow(xstar)) {
  tmp1 <- get.xprime1(xstar[i,])
  tmp2 <- get.xprime2(xstar[i,])
  xprime1 <- rbind(xprime1, tmp1)
  xprime2 <- rbind(xprime2, tmp2)
  yprime.init1 <- c(yprime.init1, yprime1(tmp1[, 1], tmp1[, 2]))
  yprime.init2 <- c(yprime.init2, yprime2(tmp2[, 1], tmp2[, 2]))
}
points(x = xprime1[,1], y = xprime1[,2], col = "red", pch = 16)
points(x = xprime2[,1], y = xprime2[,2], col = "red", pch = 16)

given <- list(x = x,
              xstar = xstar,
              xprime = rbind(xprime1, xprime2),
              y = y - mean(y))

initial.values <- list(l = c(0.5, 0.5), sig2 = 7,
                       ystar = ytrue(xstar[,1], xstar[,2]) - mean(y) +
                         rnorm(nrow(xstar), 0, 0.5),
                       yprime = c(yprime.init1, yprime.init2))

true.values <- list(ystar = ytrue(xstar[,1], xstar[,2]) - mean(y),
                    yprime = c(yprime.init1, yprime.init2))

rm(x, y, tmp1, tmp2, xprime1, xprime2, xstar,
   yprime.init1, yprime.init2, i, get.xprime1, get.xprime2)
