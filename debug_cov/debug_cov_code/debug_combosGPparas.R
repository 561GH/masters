#test.l <- abs(rmvn(6, m = c(0, 0), S = 0.1*diag(2))) 
test.l <- as.matrix(expand.grid(seq(0.02, 1, length.out = 30), seq(0.02, 1, length.out = 30)))
# sig2 values don't matter because of determinant properties
test.sig2 <-1#rchisq(6, 5)

res <- matrix(NA, nrow = nrow(test.l), ncol = length(test.sig2))  # column when length(sig2) = 1
dets <- matrix(NA, nrow = nrow(test.l), ncol = 19)  # store determinants; need length(sig2) = 1
svds <- matrix(NA, nrow = nrow(test.l), ncol = 18)  # store eigenvalues
colnames(dets) <- paste("det", 1:5, sep = "")
attach(given)
for (i in 1:nrow(test.l)) {
  for (j in 1:length(test.sig2)) {
    l <- test.l[i,]
    sig2 <- test.sig2[j]
    
    nugget1 <- nugget2 <- 0
    
    R <- covMatrix(X1 = x, X2 = x, sig2 = sig2, l = l,
                   covar.fun = "matern") +
      diag(nugget1, nrow(x))
    # Linv <- solve( t(chol(R)) )  # recall need transpose to match std chol def
    # Rinv <- t(Linv) %*% Linv
    
    n1 <- nrow(xprime) / 2  # number of derivs. in d1
    xprime1 <- xprime[1:n1,]
    xprime2 <- xprime[-(1:n1),]
    
    S11 <- covMatrix(X1 = xstar, X2 = xstar,
                     sig2 = sig2, l = l, covar.fun = "matern")
    S22 <- covMatrix(X1 = xprime1, X2 = xprime1,
                     sig2 = sig2, l = l, covar.fun = "matern2", d1 = 1, d2 = 1)
    S33 <- covMatrix(X1 = xprime2, X2 = xprime2,
                     sig2 = sig2, l = l, covar.fun = "matern2", d1 = 2, d2 = 2)
    
    S21 <- covMatrix(X1 = xprime1, X2 = xstar, # CAREFUL SEE PAPER FOR ARG. ORDER
                     sig2 = sig2, l = l, covar.fun = "matern1", d1 = 1)
    S31 <- covMatrix(X1 = xprime2, X2 = xstar, # CAREFUL SEE PAPER FOR ARG. ORDER
                     sig2 = sig2, l = l, covar.fun = "matern1", d1 = 2)
    
    S32 <- covMatrix(X1 = xprime2, X2 = xprime1,
                     sig2 = sig2, l = l, covar.fun = "matern2", d1 = 2, d2 = 1)
    
    R.xstarxprime <- rbind(cbind(S11, t(S21), t(S31)),
                           cbind(S21,   S22,  t(S32)),
                           cbind(S31,   S32,    S33)) +
      diag(rep(nugget2, nrow(xstar) + nrow(xprime)))
    
    # check determinants
    dets[i,1] <- det(S11)
    dets[i,2] <- det(S22)
    dets[i,3] <- det(S33)
    dets[i,4] <- det(R.xstarxprime[1:2, 1:2])
    dets[i,5] <- det(R.xstarxprime)
    dets[i,6] <- det(S21 %*% t(S21))
    dets[i,7] <- det(S31 %*% t(S31))
    dets[i,8] <- det(S32 %*% t(S32))
    
    dets[i,9] <- det(t(S21) %*% (S21))
    dets[i,10] <- det(t(S31) %*% (S31))
    dets[i,11] <- det(t(S32) %*% (S32))
    
    # check eigenvalues via svd
    svds[i,1] <- min(svd(S11)$d)
    svds[i,2] <- min(svd(S22)$d)
    svds[i,3] <- min(svd(S33)$d)
    svds[i,5] <- min(svd(R.xstarxprime)$d)
    svds[i,6] <- min(svd(S21)$d)
    svds[i,7] <- min(svd(S31)$d)
    svds[i,8] <- min(svd(S32)$d)
    
    svds[i,9] <- min(svd(t(S21))$d)
    svds[i,10] <- min(svd(t(S31))$d)
    svds[i,11] <- min(svd(t(S32))$d)
    
    # CAREFUL SEE PAPER FOR ARG. ORDER
    r0 <- covMatrix(X1 = xstar, X2 = x,
                    sig2 = sig2, l = l,
                    covar.fun = "matern")
    r1 <- covMatrix(X1 = xprime1, X2 = x,
                    sig2 = sig2, l = l,
                    covar.fun = "matern1", d1 = 1)
    r2 <- covMatrix(X1 = xprime2, X2 = x,
                    sig2 = sig2, l = l,
                    covar.fun = "matern1", d1 = 2)
    r.xstarprime.x <- rbind(r0,
                            r1, 
                            r2)
    
    dets[i,12] <- det(t(r0) %*% (r0))
    dets[i,13] <- det((r0) %*% t(r0))
    dets[i,14] <- det(t(r1) %*% (r1))
    dets[i,15] <- det((r1) %*% t(r1))
    dets[i,16] <- det(t(r2) %*% (r2))
    dets[i,17] <- det((r2) %*% t(r2))
    
    dets[i,18] <- det(t(r.xstarprime.x) %*% (r.xstarprime.x))
    dets[i,19] <- det((r.xstarprime.x) %*% t(r.xstarprime.x))
    
    svds[i,12] <- min(svd(t(r0))$d)
    svds[i,13] <- min(svd((r0) )$d)
    svds[i,14] <- min(svd(t(r1))$d)
    svds[i,15] <- min(svd((r1) )$d)
    svds[i,16] <- min(svd(t(r2))$d)
    svds[i,17] <- min(svd((r2) )$d)
    
    svds[i,18] <- min(svd(r.xstarprime.x)$d)
    
    m <- r.xstarprime.x %*% solve(R, y) #Rinv %*% y
    S <- R.xstarxprime -
      r.xstarprime.x %*% solve(R, t(r.xstarprime.x)) + #Rinv %*% t(r.xstarprime.x) +
      diag(rep(0, nrow(xstar) + nrow(xprime)))
    S <- ( S + t(S) ) / 2
    
    res[i, j] <- min(eigen(S)$values)
    
  }
}
detach(given)

neg <- which(res < 0)
plot(x = test.l[neg,1], y = test.l[neg,2], col = "red", pch = 16, xlab = "l1", ylab = "l2")
points(x = test.l[-neg,1], y = test.l[-neg,2], col = "black", pch = 16)
abline(v = 0.3, col = "blue")
abline(h = 0.3, col = "blue")

# Consider determinants
zz <- c(-1e-6, 1e-6)
par(mfrow = c(4,5))
apply(dets, MARGIN = 2, plot, ylim = zz, pch = 16, col = "red")
par(mfrow = c(1,1))

# compare S22
# all.equal(S22, matsecder(x = xprime1, y = xprime1, l = l, h = 1, k = 1))
zz <- c(-1e-3, 1e-3)
par(mfrow = c(4,5))
apply(svds, MARGIN = 2, plot, ylim = zz, pch = 16, col = "red")
par(mfrow = c(1,1))
