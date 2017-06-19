# Numerical instability in covariance matrix

dmvn <- function(y, m, S, log = FALSE) {

  y <- as.numeric(y)  # make it not a column
  m <- as.numeric(m)
  n <- length(y)

  # Error handling (from mvnorm)
  if (!missing(m)) {
    if (!is.null(dim(m)))
      dim(m) <- NULL
    if (length(m) != n)
      stop("m and S have non-conforming size")
  }
  if (!missing(S)) {
    if (n != ncol(S))
      stop("y and S have non-conforming size")
    if (!isSymmetric(S, tol = sqrt(.Machine$double.eps),
                     check.attributes = FALSE))
      stop("S must be a symmetric matrix")
  }

  y <- cbind(y - m)  # center y, make it a column

  L <- t( chol(S) )  # lower triangular
  alpha <- cbind( backsolve( t(L), forwardsolve(L, y)) )

  loglik <- -0.5 * t(y) %*% alpha - sum(log(diag(L))) - n/2 * log(2*pi)
  loglik <- as.numeric(loglik)

  if (log) {
    return(loglik)
  } else {
    return(exp(loglik))
  }
}

rmvn <- function(N, m, S) {

  m <- as.numeric(m)
  n <- length(m)

  # Error handling (from mvnorm)
  if (!missing(S)) {
    if (n != ncol(S))
      stop("S has non-conforming size")
    if (!isSymmetric(S, tol = sqrt(.Machine$double.eps),
                     check.attributes = FALSE))
      stop("S must be a symmetric matrix")
  }

  z <- matrix(rnorm(n * N), nrow = N, ncol = n)  # N observations = num. rows
  L <- t( chol(S) )

  trans <- function(vec) {
    matrix(L %*% vec + m, nrow = 1, ncol = n)
  }

  return( t(apply(z, MARGIN = 1, FUN = trans)) )
}
