% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulations.R
\name{simulate_tau}
\alias{simulate_tau}
\title{Simulate covariance matrix from a tau-equivalent model.}
\usage{
simulate_tau(
  n,
  k,
  lambda = 1,
  sigma = rep(1, k),
  latent = stats::rnorm,
  error = stats::rnorm
)
}
\arguments{
\item{n}{Number of observations.}

\item{k}{Size of matrix or number of testlets.}

\item{lambda}{Scalar factor loading.}

\item{sigma}{Vector of error standard deviations.}

\item{latent}{Distribution of the latent variable.}

\item{error}{Distribution of the error variable.}
}
\value{
Covariance matrix.
}
\description{
Simulate covariance matrix from a tau-equivalent model.
}