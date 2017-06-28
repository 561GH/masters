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

  # xi != xj ###################################################################
  g.xixj <- (1 + ltrans * abs(xi - xj) + 1/3 * (ltrans * abs(xi - xj))^2) *
    exp(-ltrans * abs(xi - xj))

  return( prod( g.xixj ) )

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

  # xi != xj ###################################################################
  ## for lambda = 5/2
  g.xixj <- 1/3 * ltrans^2 * abs(xi - xj) *
    (1 + ltrans * abs(xi - xj)) *
    exp( -ltrans * abs(xi - xj) ) *
    sign(xi - xj)

  return( prod(g.xixj) )

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

  # xi != xj ###################################################################
  ## for lambda = 5/2
  g.xixj <- 1/3 * ltrans^2 * (1 + ltrans * abs(xi - xj) -
                                0.5 * (ltrans * abs(xi - xj))^2 ) *
    exp( -ltrans * abs(xi - xj) ) * sign(xi - xj)

  return( prod(g.xixj) )

}
