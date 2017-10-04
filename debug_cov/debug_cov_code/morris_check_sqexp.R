C12 <- C[1:3,4:6]

Rx <- function(x, l) {
  exp(- x^2 / l)
}

Rxd <- function(x, l) {
  -2 * x * exp(- x^2 / l) / l
}

Rxdd <- function(x, l) {
  (-2 / l  + 4 * x^2 / l^2) * exp(- x^2 / l)  # morris
  #(2 / l  - 4 * x^2 / l^2) * exp(- x^2 / l)  # mine
}

t1 <- x[1,]
t2 <- x[2,]
t3 <- x[3,]
C12.13 <- Rx(0.268, l = l[1]) * Rxd(1, l = l[2])
C12.13
C12.22 <- Rxdd(0.268, l = l[1]) * Rx(1, l = l[2])
C12.22

covMatrix(X1 = rbind(t1), X2 = rbind(t2), l = l, sig2 = 1, 
          covar.fun = "sqexp2", d1 = 1, d2 = 1)  # mine is different by a sign
