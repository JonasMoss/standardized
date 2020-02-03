#' Variants of coefficient alpha
#'
#' The ordinary, standardized and sigma coefficient alpha.
#'
#' @param S A positive definite covariance matrix.
#' @return Cofficient alpha, standardized, orsigma coefficient alpha.
#' @name alpha
#' @export

alpha <- function(S) {
  k <- nrow(S)
  k / (k - 1) * (1 - tr(S) / sum(S))
}

#' @rdname alpha
#' @export
alpha_std <- function(S) {
  R <- stats::cov2cor(S)
  k <- nrow(S)
  k / (k - 1) * (1 - k / sum(R))
}

#' @rdname alpha
#' @export
alpha_sigma <- function(S) {
  stop("Function not implemented yet.")
  R <- stats::cov2cor(S)
  k <- nrow(S)
  k / (k - 1) * (1 - k / sum(R))
}

#' Variants of the omega coefficient
#'
#' The congeneric reliability, coefficient H, standardized reliability, and
#'    sigma reliability.
#'
#' @param lambda Vector of loadings.
#' @param sigma Vector of standard deviations.
#' @return The congeneric reliability, coefficient H, standardized reliability,
#'    or sigma reliability.
#' @name omega
#' @export
omega <- function(lambda, sigma) {
  k <- length(lambda)
  a <- k * mean(abs(lambda))^2 / mean(sigma^2)
  a / (a + 1)
}

#' @rdname omega
#' @export
omega_h <- function(lambda, sigma) {
  a <- sum(lambda^2 / sigma^2)
  a / (a + 1)
}

#' @rdname omega
#' @export
omega_sigma <- function(lambda, sigma) {
  k <- length(lambda)
  a <- k * mean(abs(lambda) / sigma)^2
  a / (a + 1)
}

#' @rdname omega
#' @export
omega_std <- function(lambda, sigma) {
  k <- length(lambda)
  a <- k * mean(abs(lambda) / sqrt(lambda^2 + sigma^2))^2 /
    mean(sigma^2 / (lambda^2 + sigma^2))
  a / (a + 1)
}