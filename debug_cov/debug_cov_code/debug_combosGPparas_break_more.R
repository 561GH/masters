test.l <- as.matrix(expand.grid(seq(0.02, 1, length.out = 30), seq(0.02, 1, length.out = 30)))
# sig2 values don't matter because of determinant properties
test.sig2 <-1

res1 <- matrix(NA, nrow = nrow(test.l), ncol = length(test.sig2))  # column when length(sig2) = 1
res2 <- res3 <- res4 <- res1

attach(given)

for (i in 1:nrow(test.l)) {
  for (j in 1:length(test.sig2)) {
    l <- test.l[i,]
    sig2 <- test.sig2[j]
    
    n1 <- nrow(xprime) / 2  # number of derivs. in d1
    xprime1 <- xprime[1:n1,]
    xprime2 <- xprime[-(1:n1),]
    
    # PART ONE #################################################################
    ## K is 4 x 4 blocks
    k11 <- covMatrix(X1 = xstar, X2 = xstar, sig2 = sig2, l = l,
                     covar.fun = "matern")
    k44 <- covMatrix(X1 = x, X2 = x, sig2 = sig2, l = l,
                     covar.fun = "matern")
    k41 <- covMatrix(X1 = x, X2 = xstar, sig2 = sig2, l = l,
                     covar.fun = "matern")
    
    k22 <- covMatrix(X1 = xprime1, X2 = xprime1, sig2 = sig2, l = l,
                     covar.fun = "matern2", d1 = 1, d2 = 1)
    k33 <- covMatrix(X1 = xprime2, X2 = xprime2, sig2 = sig2, l = l,
                     covar.fun = "matern2", d1 = 2, d2 = 2)
    k32 <- covMatrix(X1 = xprime2, X2 = xprime1, sig2 = sig2, l = l,
                     covar.fun = "matern2", d1 = 2, d2 = 1)
    

    k21 <- covMatrix(X1 = xprime1, X2 = xstar, sig2 = sig2, l = l,
                     covar.fun = "matern1", d1 = 1)
    k31 <- covMatrix(X1 = xprime2, X2 = xstar, sig2 = sig2, l = l,
                     covar.fun = "matern1", d1 = 2)
    k42 <- covMatrix(X1 = x, X2 = xprime1, sig2 = sig2, l = l,
                     covar.fun = "matern1", d2 = 1)
    k43 <- covMatrix(X1 = x, X2 = xprime2, sig2 = sig2, l = l,
                     covar.fun = "matern1", d2 = 2)
    
    K <- rbind( cbind(k11, t(k21), t(k31), t(k41)),
                cbind(k21,   k22 , t(k32), t(k42)),
                cbind(k31,   k32 ,   k33 , t(k43)),
                cbind(k41,   k42 ,   k43 ,   k44) )
    
    S11 <- k11 #K[1, 1]
    S12 <- cbind(t(k21), t(k31), t(k41)) #K[1, 2:4]
    S22 <- rbind( cbind(k22 , t(k32), t(k42)),
                  cbind(k32 ,   k33 , t(k43)),
                  cbind(k42 ,   k43 ,   k44) ) 
    S22 <- S22 #+ diag(rep(1e-4, nrow(S22)))
    
    # Linv <- solve( t(chol(S22)) )  # recall need transpose to match std chol def
    # S22inv <- t(Linv) %*% Linv
    # Sig1 <- S11 - S12 %*% S22inv %*% t(S12) 
    Sig1 <- S11 - S12 %*% solve(S22, t(S12)) 
    Sig1 <- (Sig1 + t(Sig1)) / 2
    
    #mu1 <- S12 %*% solve(S22, y)
    
    res1[i, j] <- min(eigen(Sig1)$values)
    
    # PART TWO #################################################################
    ## kick out the 1st row and col. of K
    # K <- rbind( cbind(k11, t(k21), t(k31), t(k41)),
    #             cbind(k21,   k22 , t(k32), t(k42)),
    #             cbind(k31,   k32 ,   k33 , t(k43)),
    #             cbind(k41,   k42 ,   k43 ,   k44) )
    S11 <- rbind( cbind(k22 , t(k32)),
                  cbind(k32 ,   k33) )
    S12 <- rbind( t(k42),
                  t(k43) )
    S22 <- k44
    
    # Linv <- solve( t(chol(S22)) )  # recall need transpose to match std chol def
    # S22inv <- t(Linv) %*% Linv
    # Sig2 <- S11 - S12 %*% S22inv %*% t(S12)
    Sig2 <- S11 - S12 %*% solve(S22, t(S12)) 
    Sig2 <- (Sig2 + t(Sig2)) / 2
    
    #mu1 <- S12 %*% solve(S22, y)
    
    res2[i, j] <- min(eigen(Sig2)$values)
    
    # PART THREE ###############################################################
    # K <- rbind(cbind(  k44 ,   k42 ,   k43 , k41),
    #            cbind(t(k42),   k22 , t(k32), k21),
    #            cbind(t(k43),   k32 ,   k33 , k31),
    #            cbind(t(k41), t(k21), t(k31), k11) )
    S11 <- k44 
    S12 <- cbind(k42 ,   k43 , k41) 
    S22 <- rbind( cbind(  k22 , t(k32), k21),
                  cbind(  k32 ,   k33 , k31),
                  cbind(t(k21), t(k31), k11) ) 
    S22 <- S22 #+ diag(rep(1e-4, nrow(S22)))
    
    # Linv <- solve( t(chol(S22)) )  # recall need transpose to match std chol def
    # S22inv <- t(Linv) %*% Linv
    # Sig1 <- S11 - S12 %*% S22inv %*% t(S12) 
    Sig3 <- S11 - S12 %*% solve(S22, t(S12)) 
    Sig3 <- (Sig3 + t(Sig3)) / 2
    
    #mu1 <- S12 %*% solve(S22, y)
    
    res3[i, j] <- min(eigen(Sig3)$values)
    
    
    # PART FOUR ################################################################

    K4 <- rbind( cbind(k11, t(k21), t(k31)),
                 cbind(k21,   k22 , t(k32)),
                 cbind(k31,   k32 ,   k33 ))
  
    res4[i, j] <- min(eigen(K4)$values)
    
  }
}
detach(given)

par(mfrow = c(2, 2))
neg1 <- which(res1 < 0)
plot(x = test.l[neg1,1], y = test.l[neg1,2], col = "red", pch = 16, cex = 2,
     xlab = "l1", ylab = "l2", xlim = c(0, 1), ylim = c(0, 1),
     main = "res1")
points(x = test.l[-neg1,1], y = test.l[-neg1,2], col = "black", pch = 16, cex = 2)

neg2 <- which(res2 < 0)
plot(x = test.l[neg2,1], y = test.l[neg2,2], col = "red", pch = 18, cex = 2,
       xlab = "l1", ylab = "l2", xlim = c(0, 1), ylim = c(0, 1),
     main = "res2")
points(x = test.l[-neg2,1], y = test.l[-neg2,2], col = "black", pch = 18, cex = 2)
#par(mfrow = c(1, 1))

#par(mfrow = c(1, 2))
neg3 <- which(res3 < 0)
plot(x = test.l[neg3,1], y = test.l[neg3,2], col = "red", pch = 16, cex = 2,
     xlab = "l1", ylab = "l2", xlim = c(0, 1), ylim = c(0, 1),
     main = "res3")
points(x = test.l[-neg3,1], y = test.l[-neg3,2], col = "black", pch = 16, cex = 2)

neg4 <- which(res4 < 0)
plot(x = test.l[neg4,1], y = test.l[neg4,2], col = "red", pch = 18, cex = 2,
     xlab = "l1", ylab = "l2", xlim = c(0, 1), ylim = c(0, 1),
     main = "res4")
points(x = test.l[-neg4,1], y = test.l[-neg4,2], col = "black", pch = 18, cex = 2)
par(mfrow = c(1, 1))


