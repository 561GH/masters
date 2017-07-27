x <- cbind(rnorm(5))


matcor(x = x, y = x, l = l)
matcorder(x = x, y = x, l = l, h = 1)
matsecder(x = x, y = x, l = l, h = 1, k = 1)

theta <- sqrt(5) / l
theta^2 / 3

matcorder(x = x, y = cbind(1:4), l = l, h = 1)
covMatrix(X = x, X2 = cbind(1:4), l = l, sig2 = s, covar.fun = r.matern1)

matcorder(x = x, y = x, l = l, h = 1)
covMatrix(X = x, X2 = x, l = l, sig2 = s, covar.fun = r.matern1)

###
X1 <- matrix(1:5, nrow = 5)
X2 <- matrix(0.5*1:3, nrow = 3)
outer(as.numeric(X1), as.numeric(X2), "-")  # X1 - X2

test <- covMatrix(X1 = X1, X2 = X2, sig2 = 1, l = 1, covar.fun = "matern")
test.sg <- matcor(x = X1, y = X2, l = 1)
all.equal(test, test.sg)  # same

test <- covMatrix(X1 = X1, X2 = X2, sig2 = 1, l = 1,
                  covar.fun = "matern1", d1 = 1)
test.sg <- matcorder(x = X1, y = X2, l = 1, h = 1)  # I think hers is wrong
all.equal(test, test.sg)

test <- covMatrix(X1 = X1, X2 = X2, sig2 = 1, l = 1,
                  covar.fun = "matern2", d1 = 1, d2 = 1)
test.sg <- matsecder(x = X1, y = X2, l = 1, h = 1, k = 1)
all.equal(test, test.sg)  # same

###
X1 <- matrix(1:10, nrow = 5)
X2 <- matrix(0.5*1:6, nrow = 3)
l <- c(1, 1)
outer(X1, X2, "-")  # too hard to interpret... use for loop

test <- covMatrix(X1 = X1, X2 = X2, sig2 = 1, l = l, covar.fun = "matern")
test.sg <- matcor(x = X1, y = X2, l = l)
all.equal(test, test.sg)  # same

test <- covMatrix(X1 = X1, X2 = X2, sig2 = 1, l = l,
                  covar.fun = "matern1", d1 = 1)
test.sg <- matcorder(x = X1, y = X2, l = l, h = 1)  # I think hers is wrong
all.equal(test, test.sg)

test <- covMatrix(X1 = X1, X2 = X2, sig2 = 1, l = l,
                  covar.fun = "matern2", d1 = 2, d2 = 2)
test.sg <- matsecder(x = X1, y = X2, l = l, h = 2, k = 2)
all.equal(test, test.sg)  # same

test <- covMatrix(X1 = X1, X2 = X2, sig2 = 1, l = l,
                  covar.fun = "matern2", d1 = 2, d2 = 1)
test.sg <- matsecder(x = X1, y = X2, l = l, h = 2, k = 1)  # I think hers is wrong
all.equal(test, test.sg)

### test symmetry
test <- covMatrix(X1 = X1, X2 = X1, sig2 = 1, l = l, covar.fun = "matern")
test.sg <- matcor(x = X1, y = X1, l = l)
all.equal(test, test.sg)

test <- covMatrix(X1 = X1, X2 = X1, sig2 = 1, l = 1, covar.fun = "matern1", d1 = 1)
test.sg <- matcorder(x = X1, y = X1, l = 1, h = 1)  # I think hers is wrong
all.equal(test, test.sg)
