% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MCMC_unconstrained.R
\name{MCMC_unconstrained}
\alias{MCMC_unconstrained}
\title{MCMC sampling for GP with unconstrained derivatives}
\usage{
MCMC_unconstrained(x, y_obs, xd, newx, N = 5000, burn = 3000,
  initial_values)
}
\arguments{
\item{x}{Evaluation sites}

\item{y_obs}{Evaluations}

\item{xd}{virtual derivative sites}

\item{newx}{prediction sites}

\item{N}{Number of MCMC iterations}

\item{burn}{burnin sample size}

\item{initial_values}{a named list of initial values for the covariance parameters \code{l} 
and \code{sig2}}
}
\value{
list of revived samples for each parameter
}
\description{
MCMC sampling for GP with unconstrained derivatives
}
