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
  k <- nrow(S)
  #stop("Not implemented")
  lambda_sq = (sum(S - diag(diag(S))))/(k * (k - 1))
  sigma_sq = diag(S) - lambda_sq

  xi <- (diag(1/sigma_sq)^(1/2)) %*% S %*% (diag(1/sigma_sq)^(1/2))
  k <- nrow(S)
  k / (k - 1) * (1 - tr(xi) / sum(xi))
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

#' The bias term for the weighted alpha.
#'
#' The bias term is described in the paper.
#'
#' @param lambda Vector of loadings.
#' @param sigma Vector of standard deviations.
#' @param w Vector of weights.
#' @return The bias term, i.e. B such that omega = alpha + B.

bias = function(lambda, sigma, w = rep(1, length(lambda))) {
  k = length(lambda)
  k / (k - 1) * (- crossprod(w, lambda)^2/k + crossprod(w^2, lambda^2)) /
    (crossprod(w, lambda)^2 + crossprod(w^2, sigma^2))
}
