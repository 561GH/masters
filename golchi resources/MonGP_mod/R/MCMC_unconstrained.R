#' MCMC sampling for GP with unconstrained derivatives
#' 
#'
#' @param x Evaluation sites
#' @param y_obs Evaluations
#' @param xd virtual derivative sites
#' @param newx prediction sites
#' @param N Number of MCMC iterations
#' @param burn burnin sample size
#' @param initial_values a named list of initial values for the covariance parameters \code{l} 
#' and \code{sig2}
#' @export
#' @return list of revived samples for each parameter 
#' 


MCMC_unconstrained = function(x, y_obs, xd, newx, N = 5000, burn = 3000, initial_values){
  M1 <- 2 
  nu <- Inf 
  l <- c()
  sig2 <- c()
  l[1] <- initial_values$l
  sig2[1] <- initial_values$sig2
  n = nrow(newx)
  m = nrow(xd)
  y <- array(dim = c(N,n))
  yd <- array(dim = c(N,m))
  R <- matcor(x, x, l[1])
  s <- solve(t(chol(R)))
  Rinv <- t(s) %*% s
  c0 <- matcor(x, newx, l[1]) 
  c1 <- matcorder(x, xd, l[1], 1)
  C <- cbind(c0, c1)
  c00 <- matcor(newx, newx, l[1]) 
  c01 <- matcorder(newx, xd, l[1], 1)
  c11 <- matsecder(xd, xd, l[1], 1, 1)
  M <- rbind(cbind(c00, c01), cbind(t(c01), c11))
  covmat <- sig2[1] * (M - t(C) %*% Rinv %*% C) 
  Rd <- covmat[(n + 1):(n + m), (n + 1):(n + m)]
  Rdinv <- ginv(Rd)
  Rdq <- 0.001 * Rd
  cc <- t(as.matrix(covmat[1:n, (n + 1):(n + m)]))
  RR <- covmat[1:n, 1:n]
  mean <- t(C) %*% Rinv %*% y_obs
  mean1 <- mean[(n + 1):(n + m)]
  mean1[mean1 < 0] <- 0
  yd[1,] <- rmnorm(1, mean1, Rdq)
  y[1, ] <- rmnorm(1, mean[1:n], RR)
  a <- 0
  b <- 0
  hh <- 0
  counter <- 0
  lpden <- c()
  lshare <- lpart(l[1], x, y_obs, newx)
  corrmat <- lshare$corrmat
  mean <- lshare$mvec
  R <- lshare$R
  covmat <- sig2 * corrmat
  lpden[1] <- lpost(l[1], sig2[1], y[1,], yd[1,], corrmat, mean, R, nu = Inf, y_obs)
  eps <- 1e-5
  for (j in 2:N) {
    
    l[j] <- l[j - 1]
    sig2[j] <- sig2[j - 1]
    y[j, ] <- y[j - 1, ]
    yd[j,] <- yd[j-1,]
    lpden[j] <- lpden[j-1]
    #############
    repeat { 
      delta <- rnorm(1, 0, 0.07)
      newl <- l[j] + delta
      if (newl> .05) 
        break
    }
    #the following vector/matrices are reused for evaluating the posterior at each iteration 
    newlshare <- lpart(newl, x, y_obs, newx)
    newcorrmat <- (t(newlshare$corrmat) + newlshare$corrmat)/2 + diag(rep(eps, (n + m)))
    newmean <- newlshare$mvec
    newR <- newlshare$R
    newcovmat <- sig2[j] * newcorrmat
    newcovmat <- (t(newcovmat) + newcovmat)/2 + diag(rep(eps, (n + m)))
    lpnum <- lpost(newl, sig2[j], y[j,], yd[j,], newcorrmat, newmean, newR, nu = Inf, y_obs)
    ratio <- (sum(log(pnorm(l[j], 0, 0.05))) + lpnum) - (sum(log(pnorm(newl, 0, 0.05))) + lpden[j])
    p <- min(1, exp(ratio))
    r <- runif(1)
    if (r <= p) {
      l[j] <- newl
      a <- a + 1
      lpden[j] <- lpnum
      corrmat<-newcorrmat
      mean<-newmean
      R<-newR
    }
    ################
    repeat {
      newsig2 <- rchisq(1, sig2[j])
      if (newsig2 > 1.7) 
        break
    }
    lpnum <- lpost(l[j], newsig2, y[j,], yd[j,], corrmat, mean, R, nu = Inf, y_obs)
    ratio <- (log(dchisq(sig2[j], newsig2)) + lpnum) - (log(dchisq(newsig2, sig2[j])) + lpden[j])
    p <- min(1, exp(ratio))
    r <- runif(1)
    if (r <= p) {
      sig2[j] <- newsig2
      b <- b + 1
      lpden[j] <- lpnum
    }
    ##################
    #generate y and y'
    covmat <- sig2[j] * corrmat
    covmat <- (t(covmat) + covmat)/2 + diag(rep(eps, n + m))
    new <- rmnorm(1, mean, covmat)
    y[j,] <- new[1:n]
    yd[j,] <- new[(n+1):(n+m)]
    lpden[j] <- lpost(l[j], sig2[j], y[j,], yd[j,], corrmat, mean, R, nu = Inf, y_obs)
    #################################################################
  }
  out <- list(l = l[(burn + 1):N], sig2 = sig2[(burn + 1):N], y = y[(burn + 1):N,], 
              yd = yd[(burn + 1):N,], lpden = lpden[(burn + 1):N])
  return(out)
}

