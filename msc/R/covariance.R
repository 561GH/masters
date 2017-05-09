# R matrix for a vector of inputs
Rmat <- function(x1, x2, sig2 = 1, l = 1) {

  # r(., .) function (p6)
  r <- function(xi, xj, sig2 = sig2, ltrans = sig2) {

    # g(.): from Matern classs of covariance functions
    ## only need g = g_k because dim(x) = 1
    ## modified bessel function of the second kind
    ## https://stat.ethz.ch/pipermail/r-help/2002-December/027850.html
    g <- function(xi, xj, ltrans) {

      # 2^(1 - lambda) / gamma(lambda) *
      #   ( sqrt(2 * lambda) * abs(xi - xj) / l ) ^ lambda *
      #   besselK(sqrt(2 * lambda) * abs(xi - xj) / l, nu = lambda)

      ## version where ltrans = sqrt(2 * lambda) / l

      # FROM NOW ON l = ltrans
      lambda <- 5/2
      2^(1 - lambda) / gamma(lambda) *
        ( ltrans * abs(xi - xj) ) ^ lambda *
        besselK( ltrans * abs(xi - xj), nu = lambda)

    }

    if (xi == xj) {
      return(sig2)
    } else {
      return(sig2 * g(xi, xj, ltrans))
    }

  }

  R <- matrix(NA, nrow = length(x1), ncol = length(x2))

  for (i in 1:length(x1)) { # I know this is inefficient
    for (j in 1:length(x2)) {

      R[i, j] <- r(xi = x1[i], xj = x2[j], sig2 = sig2, l = l)

    }
  }

  return(R)
}
