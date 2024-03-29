% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/myfnct00.R
\name{matcorder}
\alias{matcorder}
\title{First derivative correlation matrix}
\usage{
matcorder(x, y, l, h)
}
\arguments{
\item{x}{Matrix whose each row is a point in \code{R^D}}

\item{y}{Matrix whose each row is a point in \code{R^D}. The number of columns must be
the same as the number of columns in \code{x}.}

\item{l}{Vector of length scale parameters for each input dimension. Length of \code{l} must
be equal to the number of columns of \code{x} and \code{y}}

\item{h}{The dimention with respect to which the derivative is taken.}
}
\value{
Correlation matrix between points in \code{x} and \code{y}.
}
\description{
\code{matcorder} returns first partial derivative of the correlation matrix between two
sets of points given by the
matrices \code{x} and \code{y} for a given length scale parameter \code{l}.
}
