################################################################################
# Covariance of yi, yj with inputs xi, xj respectively #########################
################################################################################
r <- function(xi, xj,
              sig2,
              l = NULL, ltrans = NULL, lambda = 5/2) {

  # xi, xj: d dimensional vectors
  # sig2: constant variance parameter; non-negative real number
  # l: d dimensional length scale parameters; non-negative
  # ltrans: \sqrt(2 * \lambda) / l; d-dimensional non-negative vector
  # lambda: non-negative parameter
  # RETURN: r(xi, xj) defined on p6

  # ERROR HANDLING #############################################################
  if (missing(xi) || missing(xj)) {stop("Arguments xi, xj are missing.")}
  if (length(xi) != length(xj)) {stop("Vectors xi, xj must be same length.")}

  if (missing(sig2)) {stop("Constant variance parameter sig2 missing.")}

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
  if (isTRUE(all.equal(xi, xj))) {return(sig2)}

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

    2^(1 - lambda) / gamma(lambda) *
      ( ltrans * abs(xik - xjk) ) ^ lambda *
      besselK( ltrans * abs(xi - xj), nu = lambda)

  }
  # g FUNTION END ##############################################################

  # xi != xj ###################################################################
  xixjltrans <- cbind(xi, xj, ltrans)
  covxixj <- sig2 * prod( apply(xixjltrans, MARGIN = 1, FUN = g) )

  return(covxixj)

}

################################################################################
# Covariance matrix for a vector of inputs #####################################
################################################################################

Rmat <- function(x1, x2, sig2 = 1, l = 1) {

  R <- matrix(NA, nrow = length(x1), ncol = length(x2))

  for (i in 1:length(x1)) { # I know this is inefficient
    for (j in 1:length(x2)) {

      R[i, j] <- r(xi = x1[i], xj = x2[j], sig2 = sig2, l = l)

    }
  }

  return(R)
}
