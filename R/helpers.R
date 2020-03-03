#' Thurstone weights
#'
#' @param lambda Vector of loadings.
#' @param sigma Vector of standard deviations.
#' @return The Thurstone weights.

thurstone = function(lambda, sigma) {
  c(lambda/(sigma^2 * (1 + sum (lambda^2 / sigma^2))))
}

#' Trace of matrix
#' @param A A square matrix.
#' @return Trace of the matrix.
tr <- function(A) sum(diag(A))

#' Standardize parameter vectors
#'
#' The function `standardize_lambda` standardizes `lambda` and
#'    `standardize_sigma` standardizes `sigma`.
#' @param lambda Vector of loadings.
#' @param sigma Vector of standard deviations.
#' @return Standardized vector.
#' @name standardize
standardize_lambda = function(lambda, sigma) {
  checkmate::assert_atomic_vector(lambda, any.missing = FALSE)
  checkmate::assert_atomic_vector(sigma, any.missing = FALSE)
  checkmate::assert_numeric(lambda)
  checkmate::assert_numeric(sigma)
  lambda/sqrt(lambda^2 + sigma^2)
}

#' @rdname standardize
standardize_sigma = function(lambda, sigma) {
  checkmate::assert_atomic_vector(lambda, any.missing = FALSE)
  checkmate::assert_atomic_vector(sigma, any.missing = FALSE)
  checkmate::assert_numeric(lambda)
  checkmate::assert_numeric(sigma)
  sigma/sqrt(lambda^2 + sigma^2)
}
