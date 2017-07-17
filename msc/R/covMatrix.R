# Covariance matrix for a matrix of inputs #####################################
#' @import foreach
#' @export
covMatrix <- function(X1,
                      X2,
                      sig2,
                      l,
                      covar.fun,
                      d1 = 0, d2 = 0,
                      ...) {
  # Will always be of the form sig2 * matrix
  # X1: n1 by d design matrix; each row is a D-dimensional input
  # X2: n2 by d design matrix; each row is a D-dimensional input;
  # sig2: constant variance parameter; non-negative real number
  # l: length scale
  # covar.fun: c("matern", "matern1", "matern2"); assume lambda = 5/2
  #   if one of X1, X2 is the derivative set, X1 MUST correspond to the derivative
  #   set and X2 to the non-derivative set for r.matern1 (see paper)
  # d1: dimension for which the derivative is taken for X1 inputs
  # d2: dimension for which the derivative is taken for X2 inputs
  # ...: arguments for covar.fun
  # RETURN: n1 by n2 matrix of covariances

  # ERROR HANDLING #############################################################

  if (!is.matrix(X1) | !is.matrix(X2)) {stop("Argument X must be a matrix.")}
  if (!is.numeric(X1) | !is.numeric(X2)) {stop("Argument X must contain only numeric values.")}

  if (nrow(X1) < 1 | nrow(X2) < 1) {stop("X must have at least 1 rows (inputs).")}
  if (ncol(X2) != ncol(X1)) {
    stop("X and X2 must have the same number of columns")
  }
  if (ncol(X1) != length(l)) {
    stop("Length of l must be the same as number of columns in X1, X2.")
  }

  covar.fun <- match.arg(arg = covar.fun,
                         choices = c("matern", "matern1", "matern2"),
                         several.ok = FALSE)
  if ((covar.fun == "matern1") & d1 == 0) {stop("Specify d1.")}
  if ((covar.fun == "matern2") & (d1 == 0 | d2 == 0)) {stop("Specify d1, d2.")}

  # setup ######################################################################

  d <- ncol(X1)
  n1 <- nrow(X1)
  n2 <- nrow(X2)
  out <- array(1, dim = c(n1, n2))
  theta <- sqrt(5) / l

  # k: index for dimension
  # matern #####################################################################
  if (covar.fun == "matern") {

    for (k in 1:d) {
      diff.k <- outer(X1[, k], X2[, k], "-")
      td <- theta[k] * abs( diff.k )
      out <- out * ( 1 + td + td^2 / 3 ) * exp( -td)
    }

  }

  # matern1 ####################################################################
  # derivative is in first argument!!!
  # i.e. lower triangle for eqn (3.8) EX. r11(Xdelta, X)
  if (covar.fun == "matern1") {

    for (k in 1:d) {
      diff.k <- outer(X1[, k], X2[, k], "-")
      td <- theta[k] * abs(diff.k)

      if (k == d1) {
        out <- out * (- theta[k] / 3 * sign(diff.k) * (td + td^2) * exp(-td))
      }

      if (k != d1) {
        out <- out * ( 1 + td + td^2 / 3 ) * exp( -td)
      }
    }

  }

  # matern2 ####################################################################
  if (covar.fun == "matern2") {

    if (d1 == d2) {  # need second derivative
      for (k in 1:d) {
        diff.k <- outer(X1[, k], X2[, k], "-")
        td <- theta[k] * abs( diff.k )

        if (k == d1) {
          out <- out * (theta[k]^2 * (1 + td - td^2) / 3) * exp(-td)
        }
        if (k != d1) {
          out <- out * ( 1 + td + td^2 / 3 ) * exp( -td)
        }
      }
    }

    if (d1 != d2) {  # don't actually need second derivative
      for (k in 1:d){
        diff.k <- outer(X1[, k], X2[, k], "-")
        td <- theta[k] * abs( diff.k )

        if (k == d1) {
          out <- out * (- theta[k] / 3 * sign(diff.k) * (td + td^2) * exp(-td))
        }

        if (k == d2) {  # dg/dxj = -dg/dxi; so different by a negative; might be wrong be careful
          out <- out * (theta[k] / 3 * sign(diff.k) * (td + td^2) * exp(-td))
        }

        if (k != d1 & k != d2) {
          out <- out * ( 1 + td + td^2 / 3 ) * exp( -td)
        }
      }
    }

  }

  ## make sure resulting matrix is symmetric if relevant #######################
  out <- sig2 * out

  if (n1 == n2) {  # i.e. out is a square matrix
    out <- ( out + t(out) ) / 2
  }

  return(out)

}
