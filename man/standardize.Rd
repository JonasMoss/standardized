% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers.R
\name{standardize}
\alias{standardize}
\alias{standardize_lambda}
\alias{standardize_sigma}
\title{Standardize parameter vectors}
\usage{
standardize_lambda(lambda, sigma)

standardize_sigma(lambda, sigma)
}
\arguments{
\item{lambda}{Vector of loadings.}

\item{sigma}{Vector of standard deviations.}
}
\value{
Standardized vector.
}
\description{
The function \code{standardize_lambda} standardizes \code{lambda} and
\code{standardize_sigma} standardizes \code{sigma}.
}
