################################################################################
# Covariance matrix for a matrix of inputs #####################################
################################################################################
#' @import foreach
#' @export
covMatrix <- function(X,
                      X2 = NULL,
                      sig2,
                      covar.fun,
                      computing = "sequential",
                      ...) {
  # Will always be of the form sig2 * matrix
  # X: n by d design matrix; each row is a d-dimensional input
  # X2: n2 by d design matrix; each row is a d-dimensional input;
  # sig2: constant variance parameter; non-negative real number
  # covar.fun: choice of covariance function
  # computing: whether to compute sequentially or in parallel. Input is one
  #   of \code{c("sequential", "parallel")}; default is "sequential".
  # ...: arguments for covar.fun
  # RETURN: n by n covariance matrix defined on p6;
  # if X2 is not null, will calculate covar.fun(X, X2)
  # i.e. get n by n2 covariance matrix (n > 1, n2 > 0; i.e. min 2 by 1)

  # ERROR HANDLING #############################################################

  if (missing(X)) {stop("Argument X is missing.")}
  if (!is.matrix(X)) {stop("Argument X must be a matrix.")}
  if (!is.numeric(X)) {stop("Argument X must contain only numeric values.")}
  if (nrow(X) < 2) {stop("X must have at least 2 rows (inputs).")}

  if (missing(sig2)) {stop("Constant variance parameter sig2 missing.")}

  if (missing(covar.fun)) {stop("Missing covariance function.")}
  if (class(covar.fun) != "function") {stop("covar.fun must be a function.")}

  # PARALLEL COMPUTING #########################################################
  computing <- match.arg(arg = computing,
                         choices = c("sequential", "parallel"),
                         several.ok = FALSE)

  if (computing == "sequential") { # suppress irritating warning message for seq
    foreach::registerDoSEQ()
  } else if (foreach::getDoParWorkers() == 1) {
    cat("Selected parallel computing, but only 1 execution worker
        in the currently registered doPar backend.")
  }

  ##############################################################################

  if (!is.null(X2)) {

    if (ncol(X2) != ncol(X)) {
      stop("X and X2 must have the same number of columns")
    }
    if (!is.matrix(X2)) {stop("Argument X2 must be a matrix.")}
    if (!is.numeric(X2)) {stop("Argument X2 must contain only numeric values.")}

    n <- nrow(X)
    n2 <- nrow(X2)

    ## calculate entries by row
    ## make i, j globally bound -.-
    i <- j <- NULL
    entries <- foreach(i = 1:n, .combine = "cbind", .inorder = FALSE) %:%
      foreach(j = 1:n2, .combine = "cbind", .inorder = FALSE) %dopar% {
        do.call(what = covar.fun,
                args = list(xi = X[i,], xj = X2[j,], ...))
      }

    out <- matrix(entries, byrow = TRUE, nrow = n, ncol = n2)

  } else {

    n <- nrow(X)  # total number of inputs
    out <- matrix(NA, nrow = n, ncol = n)

    ## calculate off-diagonal entries column-wise
    ## (column-wise because of upper.tri() function)
    ## make i, j globally bound -.-
    i <- j <- NULL
    entries <- foreach(j = 2:n, .combine = "cbind", .inorder = FALSE) %:%
      foreach(i = 1:(j-1), .combine = "cbind", .inorder = FALSE) %dopar% {
        do.call(what = covar.fun,
                args = list(xi = X[i,], xj = X[j,], ...))
      }

    ## fill off-diagonals
    if (length(entries) > 1) { # or n == 2
      out[upper.tri(out)] <- entries
      out <- t(out)
      out[upper.tri(out)] <- entries
    } else {
      out[1, 2] <- out[2, 1] <- entries
    }

    ## fill diagonals
    diag(out) <- 1

  }

  return(sig2 * out)

}
