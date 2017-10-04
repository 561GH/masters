# Morris simple example

# y0 = y0(t1, t8)

## Design and data
### xi = t^i = (t1, t8)
x <- as.matrix(rbind(c(0, 0),
                     c(0.268, 1),
                     c(1, 0.268)))
y0 <- c(3.0489, 71.6374, 93.1663)
y1 <- c(12.197, 185.7917, 123.6169)
y8 <- c(27.4428, 64.1853, 244.4854)
table <- cbind(x, y0, y1, y8)

y <- c(table[1,3:5], table[2,3:5], table[3, 3:5])

v <- matrix(rep(c(1, 0, 0), 3), ncol = 1)

# posterior (assume mean 0)
r <- function(t) {
  t <- as.matrix(rbind(as.numeric(t)))
  r00 <- covMatrix(X1 = t, X2 = x, sig2 = 1, l = l, covar.fun = "matern")
  r01 <- -covMatrix(X1 = t, X2 = x, sig2 = 1, l = l,
                    covar.fun = "matern1", d1 = 1)
  r08 <- -covMatrix(X1 = t, X2 = x, sig2 = 1, l = l,
                    covar.fun = "matern1", d1 = 2)
  
  tmp <- cbind(t(r00), t(r01), t(r08))
  return(matrix(c(tmp[1,], tmp[2,], tmp[3,]), ncol = 1))
  
}

post.mean <- function(t) {
  # t: row vector!!
  r.t <- r(t)
  v <- matrix(rep(c(1, 0, 0), 3), ncol = 1)
  return(mu + t(y - mu*v) %*% solve(C, b = r.t))
}

post.covar <- function(t, s) {
  # t, s: row vector!!
  r.t <- r(t)
  r.s <- r(s)
  return( sigma2 * (1 - (t(r.t) %*% solve(C, b = r.s))) )
}


test.l <- as.matrix(expand.grid(seq(0.02, 1, length.out = 30), seq(0.02, 1, length.out = 30)))
res.mu <- res.cov <- res.sigma2 <- rep(NA, nrow(test.l))
res6 <- matrix(NA, nrow = nrow(test.l), ncol = length(test.sig2))
for (test in 1:nrow(test.l)) {
  l <- test.l[test,]
  
  C <- NULL
  for (i in 1:nrow(x)) {
    Ci <- NULL
    for (j in 1:nrow(x)) {
      xi <- rbind(x[i,])
      xj <- rbind(x[j,])
      K00 <- covMatrix(X1 = xi, X2 = xj, sig2 = 1, l = l, covar.fun = "matern")
      
      K10 <- covMatrix(X1 = xi, X2 = xj, sig2 = 1, l = l, 
                       covar.fun = "matern1", d1 = 1)
      K80 <- covMatrix(X1 = xi, X2 = xj, sig2 = 1, l = l, 
                       covar.fun = "matern1", d1 = 2)
      
      K81 <- covMatrix(X1 = xi, X2 = xj, sig2 = 1, l = l, 
                       covar.fun = "matern2", d1 = 2, d2 = 1)
      K11 <- covMatrix(X1 = xi, X2 = xj, sig2 = 1, l = l, 
                       covar.fun = "matern2", d1 = 1, d2 = 1)
      K88 <- covMatrix(X1 = xi, X2 = xj, sig2 = 1, l = l, 
                       covar.fun = "matern2", d1 = 2, d2 = 2)
      
      Cij <- rbind(c(K00, t(K10), t(K80)),
                   c(K10,   K11,  t(K81)),
                   c(K80,   K81,    K88))
      
      Ci <- cbind(Ci, Cij)
    }
    C <- rbind(C,
               Ci)
  }
  
  C <- (C + t(C)) / 2  
  
  n <- 3
  k <- 2
  mu <- t(v) %*% solve(C, y) / ( t(v) %*% solve(C, v) )
  mu <- as.numeric(mu)
  sigma2 <- 1 / (n * (k + 1)) * t(y - mu * v) %*% solve(C, (y - mu * v)) 
  sigma2 <- as.numeric(sigma2) 
  
  res.sigma2[test] <- sigma2
  res.mu[test] <- post.mean(t = cbind(0.5, 0.5))
  res.cov[test] <- post.covar(t = cbind(0.5, 0.5), s = cbind(0.5, 0.5))
  
  # test my problems
  t <- cbind(3, 5) / 10
  s <- cbind(5, 3) / 10
  r.t <- r(t)
  r.s <- r(s)
  # S11 <- sigma2 * (1 - (t(r.t) %*% solve(C, b = r.t)))
  # S22 <- sigma2 * (1 - (t(r.s) %*% solve(C, b = r.s)))
  # S12 <- sigma2 * (1 - (t(r.t) %*% solve(C, b = r.s)))
  
  
  S11 <- covMatrix(X1 = rbind(t, s), X2 = rbind(t, s), sig2 = 1, l = l, covar.fun = "matern")
  S12 <- rbind(t(r.t), t(r.s))
  
  S <- S11 - S12 %*% solve(C, t(S12))
  S <- (S + t(S)) / 2
  res6 <- min(eigen(S)$values)
}

# test
res <- cbind(test.l, res.cov, res.mu, res.sigma2)
res <- res[order(res.cov, decreasing = TRUE),]
res <- res[res[,3] > 0, ]
tail(res)

target.mu <- res[order(res[,4], decreasing = FALSE),]

# my problem
neg6 <- which(res6 < 0)
plot(x = test.l[neg6,1], y = test.l[neg6,2], col = "red", pch = 18, cex = 2,
     xlab = "l1", ylab = "l2", xlim = c(0, 1), ylim = c(0, 1),
     main = "res4")
points(x = test.l[-neg6,1], y = test.l[-neg6,2], col = "black", pch = 18, cex = 2)
