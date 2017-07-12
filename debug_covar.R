# Diagnosing covariance function for the derivative
## run: matcor, matcorder, matsecder
matcor = function(x, y, l){
  t = sqrt(5)/l
  M = length(t)
  n1 = nrow(x)
  n2 = nrow(y)
  out = array(1, dim = c(n1, n2))
  for (d in 1:M){
    D = outer(x[,d], y[,d], "-")
    D = t[d] * abs(D)
    out = out*(1 + D + (D ^ 2) / 3) * exp(-D)
  }
  return(out)
}

#' First derivative correlation matrix
#'
#'
#' \code{matcorder} returns first partial derivative of the correlation matrix between two
#' sets of points given by the
#' matrices \code{x} and \code{y} for a given length scale parameter \code{l}.
#' @param x Matrix whose each row is a point in \code{R^D}
#' @param y Matrix whose each row is a point in \code{R^D}. The number of columns must be
#' the same as the number of columns in \code{x}.
#' @param l Vector of length scale parameters for each input dimension. Length of \code{l} must
#' be equal to the number of columns of \code{x} and \code{y}
#' @param h The dimention with respect to which the derivative is taken.
#' @return Correlation matrix between points in \code{x} and \code{y}.
#'

matcorder = function(x, y, l, h){
  t = sqrt(5) / l
  M = length(t)
  n1 = nrow(x)
  n2 = nrow(y)
  out = array(1, dim = c(n1, n2))
  for (d in 1:M){
    D = outer(x[,d], y[,d], "-")
    I = sign(D)
    if (d == h) {
      D = t[d] * abs(D)
      out =  out * I * (t[d] * D * (1 + D) / 3) * exp(-D)}
    if (d!=h) {
      D = t[d] * abs(D)
      out = out * (1 + D + (D ^ 2) / 3) * exp(-D)}
  }
  return(out)
}

#' Second derivative correlation matrix
#'
#'
#' \code{matsecder} returns second derivative of the correlation matrix between two
#' sets of points given by the
#' matrices \code{x} and \code{y} for a given length scale parameter \code{l}.
#' @param x Matrix whose each row is a point in \code{R^D}
#' @param y Matrix whose each row is a point in \code{R^D}. The number of columns must be
#' the same as the number of columns in \code{x}.
#' @param l Vector of length scale parameters for each input dimension. Length of \code{l} must
#' be equal to the number of columns of \code{x} and \code{y}
#' @param h The dimention with respect to which the first derivative is taken.
#' @param k The dimention with respect to which the first derivative is taken.
#' @return Correlation matrix between points in \code{x} and \code{y}.
#'

matsecder = function(x, y, l, h, k){
  t = sqrt(5) / l
  M = length(t)
  n1 = nrow(x)
  n2 = nrow(y)
  out = array(1, dim = c(n1, n2))
  if (h == k) {
    for (d in 1:M){
      D = outer(x[,d], y[,d], "-")
      D = t[d] * abs(D)
      if (d == h) {
        out = out * (t[d] ^ 2 * (1 + D - D ^ 2) / 3) * exp(-D)
      }
      if (d != h) {
        out = out * (1 + D + (D ^ 2) / 3) * exp(-D)
      }
    }
  }
  if (h != k) {
    for (d in 1:M){
      D = outer(x[,d], y[,d], "-")
      I = sign(D)
      D = t[d] * abs(D)
      if (d == h | d == k) {
        out =  out * I * (t[d] * D * (-1 - D) / 3) * exp(-D)
      }
      else {
        out = out * (1 + D + (D ^ 2) / 3) * exp(-D)
      }
    }
    out = -out
  }
  return(out)}

library(mvtnorm)
ytrue <- function(x) {log(20 * x + 1)}  # x > -1/20; monotone increasing function

# x = inputs where fcn evaluated
# xstar = prediction set inputs
# xprime = derivative set inputs
# true/known model values
x <- cbind(seq(0, 1, length.out = 10))
xstar <- rbind(0.4, 0.42, 0.8, 0.82)
xprime <- rbind(0.41, 0.81)
y <- c( ytrue(x) ) - mean(ytrue(x))  # to have 0 mean

# Comparing ####################################################################
s <- 1  # sig2
l <- 0.3  # length scale

Kyy <- covMatrix(X = x, sig2 = s, covar.fun = r.matern, l = l)
Kyy.sg <- s * matcor(x = x, y = x, l = l)
all.equal(Kyy, Kyy.sg)  # yes

## first derivative
Kystarystar <- covMatrix(X = xstar, sig2 = s, covar.fun = r.matern1, l = l)
Kystarystar.sg <- s * matcorder(x = xstar, y = xstar, l = l, h = 1)
Kystarystar
Kystarystar.sg  # not symmetric?

## second derivative
Kystarystar2 <- covMatrix(X = xstar, sig2 = s, covar.fun = r.matern2, l = l)
Kystarystar2.sg <- s * matsecder(x = xstar, y = xstar, l = l, h = 1, k = 1)
Kystarystar2
Kystarystar2.sg

## whatever for now
plot(x = x, y = y, pch = 16, cex = 2)
points(x = xstar, y = ytrue(xstar) - mean(ytrue(x)), col = "red")
Linv <- solve( t(chol(Kyy)) )  # recall need transpose to match std chol def
Kyyinv <- t(Linv) %*% Linv

Kystarystar <- covMatrix(X = xstar, sig2 = s, covar.fun = r.matern, l = l)
Kystary <- covMatrix(X = xstar, X2 = x, sig2 = s, covar.fun = r.matern, l = l)

post.mean <- Kystary %*% Kyyinv %*% cbind(y)
post.covar <- Kystarystar - Kystary %*% Kyyinv %*% t(Kystary)
post.covar <- ( post.covar + t(post.covar) ) / 2
points(x = xstar, y = post.mean, col = "blue", pch = 16)

sple <- rmvnorm(3, mean = post.mean, sigma = post.covar)
for ( i in 1:nrow(sple) ) {
  points(x = xstar, y = sple[i,], col = "red", pch = 10)
}

# some sample paths
xnew <- cbind(seq(0, 1, length.out = 49))
Kynewynew <- covMatrix(X = xnew, sig2 = s, covar.fun = r.matern, l = l)
Kynewy <- covMatrix(X = xnew, X2 = x, sig2 = s, covar.fun = r.matern, l = l)

post.mean <- Kynewy %*% Kyyinv %*% cbind(y)
post.covar <- Kynewynew - Kynewy %*% Kyyinv %*% t(Kynewy)
post.covar <- ( post.covar + t(post.covar) ) / 2

sple.paths <- rmvnorm(5, mean = post.mean, sigma = post.covar)

for ( i in 1:nrow(sple.paths) ) {
  points(x = xnew, y = sple.paths[i,], type = "l", col = "blue")
}

# numerically get derivatives at xprime ########################################
yprime.true <- 20 / (20 * xprime + 1)  # too strict though
slope1 <- (sple[,2] - sple[,1]) / (xstar[2] - xstar[1])
slope2 <- (sple[,4] - sple[,3]) / (xstar[4] - xstar[3])

# now try with derivative GP
Kyprimeyprime <- covMatrix(X = xprime, sig2 = s, covar.fun = r.matern2, l = l)
Kyprimey <- covMatrix(X = xprime, X2 = x, sig2 = s, covar.fun = r.matern1, l = l)

post.mean <- Kyprimey %*% Kyyinv %*% cbind(y)
post.covar <- Kyprimeyprime - Kyprimey %*% Kyyinv %*% t(Kyprimey)
post.covar <- ( post.covar + t(post.covar) ) / 2
eigen(post.covar)$values

# mine is wrong... try shirin's
Kyprimeyprime.sg <- s * matsecder(x = xprime, y = xprime, l = l, h = 1, k = 1)
Kyprimey.sg <- s * matcorder(x = xprime, y = x, l = l, h = 1)

post.mean.sg <- Kyprimey.sg %*% Kyyinv %*% cbind(y)
post.covar.sg <- Kyprimeyprime.sg - Kyprimey.sg %*% Kyyinv %*% t(Kyprimey.sg)
post.covar.sg <- ( post.covar + t(post.covar) ) / 2
eigen(post.covar)$values


# exactly the same posterior mean = both negative values i.e. both are super wrong
# my covariance seems more wrong... eigenvalues are super negative
