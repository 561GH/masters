################################################################################
# Covariance of yi, yj with inputs xi, xj respectively #########################
################################################################################

################################################################################
# Matern #######################################################################
################################################################################
#' @export
r.matern <- function(xi, xj,
              l = NULL, ltrans = NULL, lambda = 5/2,
              ...) {

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
      ( ltransk * abs(xik - xjk) ) ^ lambda *
      besselK( ltransk * abs(xik - xjk), nu = lambda)

    ## for lambda = 5/2; used to check if the above is correct
    # (1 + ltransk * abs(xik - xjk) + 1/3 * (ltransk * abs(xik - xjk))^2) *
    #   exp(-ltransk * abs(xik - xjk))

  }
  # g FUNCTION END #############################################################

  # xi != xj ###################################################################
  xixjltrans <- cbind(xi, xj, ltrans)
  covxixj <- prod( apply(xixjltrans, MARGIN = 1, FUN = g) )

  return(covxixj)

}

################################################################################
# Matern1 ######################################################################
################################################################################
#' @export
r.matern1 <- function(xi, xj,
                     l = NULL, ltrans = NULL, lambda = 5/2,
                     ...) {

  # Only for lambda = 5/2 for now.
  #
  # xi, xj: d dimensional vectors;
  #   xi MUST correspond to the derivative set, xj MUST correspond to the
  #   non-derivative for r.matern1 (see paper)
  # l: d dimensional length scale parameters; non-negative
  # ltrans: \sqrt(2 * \lambda) / l; d-dimensional non-negative vector
  # lambda: non-negative parameter
  # RETURN: r(xi, xj) defined on p6 with derivative of matern (appendix)

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

    ## for lambda = 5/2
    1/3 * ltransk^2 * abs(xik - xjk) *
      (1 + ltransk * abs(xik - xjk)) *
      exp( -ltransk * abs(xik - xjk) ) *
      sign(xik - xjk)

  }
  # g FUNCTION END #############################################################

  # xi != xj ###################################################################
  xixjltrans <- cbind(xi, xj, ltrans)
  covxixj <- prod( apply(xixjltrans, MARGIN = 1, FUN = g) )

  return(covxixj)

}

################################################################################
# Matern2 ######################################################################
################################################################################
#' @export
r.matern2 <- function(xi, xj,
                      l = NULL, ltrans = NULL, lambda = 5/2,
                      ...) {

  # Only for lambda = 5/2 for now.
  #
  # xi, xj: d dimensional vectors
  # l: d dimensional length scale parameters; non-negative
  # ltrans: \sqrt(2 * \lambda) / l; d-dimensional non-negative vector
  # lambda: non-negative parameter
  # RETURN: r(xi, xj) defined on p6 with derivative of matern (appendix)

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

    ## for lambda = 5/2
    1/3 * ltransk^2 * (1 + ltransk * abs(xik - xjk) -
                        0.5 * (ltransk * abs(xik - xjk))^2 ) *
      exp( -ltrans * abs(xik - xjk) ) * sign(xik - xjk)

  }
  # g FUNCTION END #############################################################

  # xi != xj ###################################################################
  xixjltrans <- cbind(xi, xj, ltrans)
  covxixj <- prod( apply(xixjltrans, MARGIN = 1, FUN = g) )

  return(covxixj)

}
