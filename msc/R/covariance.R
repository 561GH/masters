################################################################################
# Covariance of yi, yj with inputs xi, xj respectively #########################
################################################################################
#' @export
r <- function(xi, xj,
              l = NULL, ltrans = NULL, lambda = 5/2) {

  # xi, xj: d dimensional vectors
  # l: d dimensional length scale parameters; non-negative
  # ltrans: \sqrt(2 * \lambda) / l; d-dimensional non-negative vector
  # lambda: non-negative parameter
  # RETURN: r(xi, xj) defined on p6

  # ERROR HANDLING #############################################################
  if (missing(xi) || missing(xj)) {stop("Arguments xi, xj are missing.")}
  if (length(xi) != length(xj)) {stop("Vectors xi, xj must be same length.")}

  if (is.null(l) & is.null(ltrans)) {
    stop("Length scale parameter l or ltrans missing.")
  }
  if (!is.null(l) & !is.null(ltrans)) {
    stop("Specify length scale parameter as either l or ltrans, not both.")
  }

  # always work with ltrans ####################################################
  if (is.null(ltrans)) {ltrans <- sqrt(2 * lambda) / l}
  if (length(xi) != length(ltrans)) {
    stop("ltrans must be same length as xi, xj.")
  }

  # if xi == xj ################################################################
  if (isTRUE(all.equal(xi, xj))) {return(1)}  # multiply by sig2 in covMatrix()

  # g FUNCTION START ###########################################################
  # g(.): from Matern classs of covariance functions
  ## modified bessel function of the second kind
  ## https://stat.ethz.ch/pipermail/r-help/2002-December/027850.html
  g <- function(input) {
    # input: 3-dimensional vector c(xik, xjk, ltransk)
    ## xik, xjk: kth entry of vectors xi, xj respectively
    ## ltransk: kth entry of ltrans as defined for r(., .); always ltrans, not l

    xik <- input[1]
    xjk <- input[2]
    ltransk <- input[3]

    ## for general lambda
    2^(1 - lambda) / gamma(lambda) *
      ( ltrans * abs(xik - xjk) ) ^ lambda *
      besselK( ltrans * abs(xik - xjk), nu = lambda)

    ## for lambda = 5/2; used to check if the above is correct
    # (1 + ltrans * abs(xik - xjk) + 1/3 * (ltrans * abs(xik - xjk))^2) *
    #   exp(-ltrans * abs(xik - xjk))

  }
  # g FUNTION END ##############################################################

  # xi != xj ###################################################################
  xixjltrans <- cbind(xi, xj, ltrans)
  covxixj <- prod( apply(xixjltrans, MARGIN = 1, FUN = g) )

  return(covxixj)

}

################################################################################
# Covariance matrix for a matrix of inputs #####################################
################################################################################
#' @import foreach
#' @export
covMatrix <- function(X,
                      sig2,
                      l = NULL, ltrans = NULL, lambda = 5/2,
                      computing = "sequential") {
  # X: n by d design matrix; each row is a d-dimensional input
  # sig2: constant variance parameter; non-negative real number
  # l: d dimensional length scale parameters; non-negative
  # ltrans: \sqrt(2 * \lambda) / l; d-dimensional non-negative vector
  # lambda: non-negative parameter
  # computing: whether to compute sequentially or in parallel. Input is one
  #   of \code{c("sequential", "parallel")}; default is "sequential".
  # RETURN: n by n covariance matrix defined on p6

  # ERROR HANDLING #############################################################

  if (missing(X)) {stop("Argument X is missing.")}
  if (!is.matrix(X)) {stop("Argument X must be a matrix.")}
  if (!is.numeric(X)) {stop("Argument X must contain only numeric values.")}
  if (nrow(X) < 2) {stop("X must have at least 2 rows (inputs).")}

  if (missing(sig2)) {stop("Constant variance parameter sig2 missing.")}

  if (is.null(l) & is.null(ltrans)) {
    stop("Length scale parameter l or ltrans missing.")
  }
  if (!is.null(l) & !is.null(ltrans)) {
    stop("Specify length scale parameter as either l or ltrans, not both.")
  }

  # always work with ltrans ####################################################
  if (is.null(ltrans)) {ltrans <- sqrt(2 * lambda) / l}
  if (ncol(X) != length(ltrans)) {
    stop("ltrans must be same dimensions as the inputs i.e. ncol(X).")
  }

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

  n <- nrow(X)  # total number of inputs
  out <- matrix(NA, nrow = n, ncol = n)

  ## calculate off-diagonal entries column-wise
  ## (column-wise because of upper.tri() function)
  ## make i, j globally bound -.-
  i <- j <- NULL
  entries <- foreach(j = 2:n, .combine = "cbind", .inorder = FALSE) %:%
    foreach(i = 1:(j-1), .combine = "cbind", .inorder = FALSE) %dopar% {
      r(X[i,], X[j,], ltrans = ltrans, lambda = lambda)
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

  return(sig2 * out)

}
