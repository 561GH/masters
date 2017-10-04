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

## prior covariance matrix: Sigma; set sigma2 = 1
### using my cov function
x <- rbind(c(0, 0),
           c(0.268, 1),
           c(1, 0.268))
l <- 1 / c(0.429, 0.467)

C <- NULL
for (i in 1:nrow(x)) {
  Ci <- NULL
  for (j in 1:nrow(x)) {
    xi <- rbind(x[i,])  # order of arguments different in morris see (1.16)
    xj <- rbind(x[j,])
    K00 <- covMatrix(X1 = xi, X2 = xj, sig2 = 1, l = l, covar.fun = "sqexp")
    
    K10 <- covMatrix(X1 = xi, X2 = xj, sig2 = 1, l = l, 
                     covar.fun = "sqexp1", d1 = 1)
    K80 <- covMatrix(X1 = xi, X2 = xj, sig2 = 1, l = l, 
                     covar.fun = "sqexp1", d1 = 2)
    
    K81 <- covMatrix(X1 = xi, X2 = xj, sig2 = 1, l = l, 
                     covar.fun = "sqexp2", d1 = 2, d2 = 1)
    K11 <- covMatrix(X1 = xi, X2 = xj, sig2 = 1, l = l, 
                     covar.fun = "sqexp2", d1 = 1, d2 = 1)
    K88 <- covMatrix(X1 = xi, X2 = xj, sig2 = 1, l = l, 
                     covar.fun = "sqexp2", d1 = 2, d2 = 2)
    
    Cij <- rbind(c(K00,  -K10,  -K80),
                 c(K10,   K11,   K81),
                 c(K80,   K81,   K88))
    
    Ci <- cbind(Ci, Cij)
  }
  C <- rbind(C,
             Ci)
}

det(C)

n <- 3
k <- 2
mu <- t(v) %*% solve(C, y) / ( t(v) %*% solve(C, v) )
mu <- as.numeric(mu)
sigma2 <- 1 / (n * (k + 1)) * t(y - mu * v) %*% solve(C, (y - mu * v))
sigma2 <- as.numeric(sigma2)

mu <- 69.15
sigma2 <- 135.47^2

# posterior (assume mean 0)
r <- function(t) {
  t <- as.matrix(rbind(as.numeric(t)))
  r00 <- covMatrix(X1 = t, X2 = x, sig2 = 1, l = l, covar.fun = "sqexp")
  r01 <- covMatrix(X1 = t, X2 = x, sig2 = 1, l = l,
                        covar.fun = "sqexp1", d2 = 1)
  r08 <- covMatrix(X1 = t, X2 = x, sig2 = 1, l = l,
                        covar.fun = "sqexp1", d2 = 2)
  
  
  tmp <- rbind(r00, r01, r08)
  return(matrix(c(tmp[,1], tmp[,2], tmp[,3]), ncol = 1) )
  
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

# test
post.mean(t = cbind(0.5, 0.5))
post.covar(t = cbind(0.5, 0.5), s = cbind(0.5, 0.5)); 2.7^2

post.mean(t = cbind(1, 1))
sqrt( post.covar(t = cbind(1, 1), s = cbind(1, 1)) )

# test contours 
t.star <- as.matrix(expand.grid(x = seq(0, 1, length.out = 10), y = seq(0, 1, length.out = 10)))
z <- apply(t.star, 1, FUN = post.mean)

library(plotly)

training <- plot_ly(x = x[, 1], y = x[, 2], z = y0, type = "contour")
training
p <- plot_ly(x = t.star[, 1], y = t.star[, 2], z = z, type = "contour")
p

# change mean to be 0 