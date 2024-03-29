% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/myfnct00.R
\name{lprioryyd}
\alias{lprioryyd}
\title{Log-prior on response and derivatives}
\usage{
lprioryyd(y, yd, mean, covmat)
}
\arguments{
\item{y}{Vector of responses.}

\item{yd}{Vector of derivative values.}

\item{mean}{Mean vector of the prior predictive distribution.}

\item{covmat}{Covariance matrix of the prior predictive distribution.}
}
\value{
Log-prior predictive density
}
\description{
evaluates log-prior predictive density for a set of responses and derivatives
}
