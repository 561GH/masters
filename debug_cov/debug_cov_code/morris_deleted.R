# res <- NULL
# for (i in 1:nrow(x)) {
#   xi <- rbind(x[i,])
#   r00 <- covMatrix(X1 = t, X2 = xi, sig2 = 1, l = l, covar.fun = "matern")
#   r01 <- -covMatrix(X1 = t, X2 = xi, sig2 = 1, l = l, 
#                     covar.fun = "matern1", d1 = 1)
#   r08 <- -covMatrix(X1 = t, X2 = xi, sig2 = 1, l = l, 
#                     covar.fun = "matern1", d1 = 2)
#   res <- c(res, r00, r01, r08)
# }
# return(res)