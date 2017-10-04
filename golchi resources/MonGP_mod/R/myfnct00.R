#' Correlation matrices
#'
#'
#' \code{matcor} returns the correlation matrix between two sets of points given by the
#' matrices \code{x} and \code{y} for a given length scale parameter \code{l}.
#' @param x Matrix whose each row is a point in \code{R^D}
#' @param y Matrix whose each row is a point in \code{R^D}. The number of columns must be
#' the same as the number of columns in \code{x}.
#' @param l Vector of length scale parameters for each input dimension. Length of \code{l} must
#' be equal to the number of columns of \code{x} and \code{y}
#' @return Correlation matrix between points in \code{x} and \code{y}.
#'

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

#' Log Multivariate Gausssian likelihood
#'
#'
#' \code{lmvnlik} returns the log pdf of the multivariate normal. There are other functions that
#' do the exact same thing. But there are some manual error preventions here that might make
#' your life easier; or not!
#' @param y Vector at which the pdf is obtained.
#' @param mean Vector mean of the multivariate Gaussian. Same length as \code{y}.
#' @param covmat covariance matrix of the multivariate Gaussian.
#' @return Log multivariate Gaussian pdf with mean \code{mean} and covariance matrix \code{covmat}
#' for vector \code{y}
#'

lmvnlik = function(y, mean, covmat){
  n = length(mean)
  s = solve(t(chol(covmat)))
  covinv = t(s) %*% s
  det = det(covinv)
  if (det < 0) det = abs(det)
  a =  -t(y - mean) %*% covinv %*% (y - mean) / 2
  b =  (n / 2) * log(2 * pi) - (0.5) * log(det)
  out =  a - b
  return(out)}

#' A Metropolis-Hastings
#'
#' Metropolis-Hastings step to draw output and derivative from the constrained GP with fixed GP parameters
#' @param initial Initial vector of outputs and derivatives.
#' @param mean Stacked vector mean of the GP and its derivative.
#' @param covmat Stacked covariance matrix.
#' @param nu Constraint parameter - controls the monotonicity.
#' @param q Scale parameter for the covariance matrix of the proposal distribution.
#' @param N Number of MCMC iterations.
#' @return A list: first element is a draw from the constrained GP. Second element is the
#' acceptance rate of the MCMC.
#'

MH = function(initial, mean, covmat, nu, q, M1, n){
  k = 0
  len = length(mean)
  z = mat.or.vec(M1, len)
  mean = c(mean)
  z[1,] = initial
  if (covmat[1,1] <= 1e-4) {
    for (i in 2:M1) {
      z[i,] = initial
      new <- ratio <- NA
      print(ratio)
    }
  } else {
    pcov = q * covmat
    for (i in 2:M1){
      z[i,] = z[i-1,]
      new = rmnorm(1, z[i,], pcov)
      g0 = sum(log(pnorm(new[(n+1):len]/nu)))
      g1 = sum(log(pnorm(z[i,][(n + 1):len] / nu)))
      ratio = lmvnlik(c(new),mean,covmat) + g0 - lmvnlik(z[i,], mean, covmat) - g1
      print(ratio)
      if (runif(1) <= min(1,exp(ratio))) {
        z[i,] = c(new)
        k = k+1
      }
    }
  }
  out = list(draw = z[M1,], acc.r = k, new = new, ratio = ratio)
  return(out)
  }

#' Log Prior on GP parameters
#'
#' Evaluates log-prior for a set of GP parameters
#' @param l Vector of length scale parameters.
#' @param sig2 Variance. scalar.
#' @param lmean The mean of the chi-squared prior on \code{l}. Default value is 1.
#' @param sigmean The mean of the chi-squared prior on \code{l}. Default value is 5.
#' @return Log-prior density at the given \code{l} and \code{sig2} pair.
#'


logprior = function(l, sig2, lmean = 1, sigmean = 5) {
  term1 = sum(log(sqrt(5) / (l ^ 2))) + sum(dchisq(sqrt(5)/l, 1, log = T))
  term2 = dchisq(sig2, 5, log = T)
  out = term1 + term2
  return(out)
}

#' Update function
#'
#' Updates the objects in sampling that are a function of the length scale \code{l}
#' @param l Vector of length scale parameters.
#' @param x Design matrix.
#' @param y Vector of responses.
#' @param newx Prediction points.
#' @return list of correlation matrix for the design points; mean vector and covariance
#' matrix of the posterior predictive
#'

lpart = function(l, x, y_obs, newx) {
  R = matcor(x, x, l)
  Rinv = ginv(R)
  c0 = matcor(x, newx, l)
  c1 = matcorder(x, xd, l, 1)
  C = cbind(c0, c1)
  c00 = matcor(newx, newx, l)
  c01 = matcorder(newx, xd, l, 1)
  c11 = matsecder(xd, xd, l, 1, 1)
  M = rbind(cbind(c00, c01), cbind(t(c01), c11))
  mean = t(C) %*% Rinv %*% y_obs
  corrmat = M - t(C) %*% Rinv %*% C
  out = list(R = R, corrmat = corrmat, mvec = mean)
  return(out)
}

#' Log-prior on response and derivatives
#'
#' evaluates log-prior predictive density for a set of responses and derivatives
#' @param y Vector of responses.
#' @param yd Vector of derivative values.
#' @param mean Mean vector of the prior predictive distribution.
#' @param covmat Covariance matrix of the prior predictive distribution.
#' @return Log-prior predictive density
#'

lprioryyd = function(y, yd, mean, covmat) {
  bigy = c(y, yd)
  out = lmvnlik(bigy, mean, covmat)
  return(out)
}

#' GP log-likelihood
#'
#' evaluates GP log-likelihood for a set of responses given variance and correlation matrix
#' @param y_obs Vector of outputs.
#' @param R Correlation matrix.
#' @param sig2 scalar variance.
#' @return log-likelihood
#'

llik = function(y_obs, R, sig2) return(lmvnlik(y_obs, rep(0, length(data)), sig2 * R))

#' Constraint function
#'
#' @param yd Vector of derivatives.
#' @param nu Constraint parameter.
#' @return Constraint value, as nu gets large returns 0 if yd is negative and 1 otherwise.
#'

cons <- function(yd, nu){
  out <- sum(log(pnorm(yd/nu)))
  return(out)
}

#' Log posterior function
#'
#' @param l Correlation parameter
#' @param sig2 Variance
#' @param y Vector of responses
#' @param yd Vector of derivatives
#' @param corrmat Posterior predictive correlation matrix
#' @param mean GP mean vector
#' @param R Likelihood correlation matrix
#' @param nu Constraint parameter
#' @return log posterior predictive value
#'

lpost <- function(l, sig2, y, yd, corrmat, mean, R, nu, y_obs) {
  eps <- 1e-5
  covmat <- sig2 * corrmat + diag(rep(eps, dim(corrmat)[[1]]))
  out <- logprior(l, sig2) + lprioryyd(y, yd, mean, covmat) + llik(y_obs, sig2, R) + cons(yd, nu)
  if (out == Inf) out <- -Inf
  return(out)
}

#' MCMC sampling within SMC
#'
#'
#' @param l Correlation parameter
#' @param sig2 Variance
#' @param y Vector of responses
#' @param yd Vector of derivatives
#' @param M Number of MCMC transition, default is 1.
#' @param nu Constraint parameter
#' @export
#' @return list of revived samples for each parameter
#'

sampling <- function(l, sig2, y, yd, M, nu) {
  ql = sd(l)
  N = length(l)
  n = ncol(y)
  m = ncol(yd)
  debug = NULL
  for (i in 1:N) {
    for (j in 1:M) {
      repeat {
        delta <- rnorm(1, 0, ql)
        newl <- l[i] + delta
        if (newl > 0)
          break
      }
      newlshare <- lpart(newl, x, y_obs, newx)
      newcorrmat <- newlshare$corrmat
      newmean <- newlshare$mvec
      newR <- newlshare$R
      newcovmat <- sig2[i] * newcorrmat
      lpnum <- lpost(newl, sig2[i], y[i,], yd[i,], newcorrmat, newmean, newR, nu, y_obs)
      lsharej <- lpart(l[i], x, y_obs, newx)
      corrmatj <- lsharej$corrmat
      meanj <- lsharej$mvec
      Rj <- lsharej$R
      covmatj <- sig2[i] * corrmatj
      lpden <- lpost(l[i], sig2[i], y[i,], yd[i,], corrmatj, meanj, Rj, nu, y_obs)
      ratio <- (log(pnorm(l[i], 0, 0.01)) + lpnum) - (log(pnorm(newl, 0, 0.01)) + lpden)
      p <- min(1, exp(ratio))
      r <- runif(1)
      if (r <= p) {
        l[i] <- newl
        #a <- a + 1
        corrmatj <- newcorrmat
        meanj <- newmean
        Rj <- newR
        covmatj <- newcovmat
        lpden <- lpnum
      }
      repeat {
        newsig2 <- rchisq(1, sig2[i])
        if (newsig2 > 0.5)
          break
      }
      lpnum <- lpost(l[i], newsig2, y[i,], yd[i,], corrmatj, meanj, Rj, nu, y_obs)
      ratio <- (log(dchisq(sig2[i], newsig2)) + lpnum) - (log(dchisq(newsig2, sig2[i]))+lpden)
      p <- min(1, exp(ratio))
      r <- runif(1)
      if (r <= p) {
        sig2[i] <- newsig2
        #b <- b + 1
      }
      covmatj <- sig2[i] * corrmatj
      covmatj <- (t(covmatj) + covmatj)/2
      sple <- MH(c(y[i,], yd[i,]), meanj, covmatj, nu, .01, M1 = 2, n)
      z <- sple$draw
      #list(draw = z[M1,], acc.r = k, new = new, ratio = ratio)
      debug <- rbind(debug, c(acc = sple$acc.r, new = sple$new, hr = sple$ratio))
      y[i,] <- z[1:n]
      yd[i,] <- z[(n+1):(n+m)]
      #k[i]<-z[[2]]
    }
  }
  out <- list(l = l, sig2 = sig2, y = y, yd = yd, debug = debug)
  return(out)
}
