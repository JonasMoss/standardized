#' Theoretical ordinal reliability
#'
#' @param lambda Vector of loadings.
#' @param sigma Vector of standard deviations.
#' @param cuts List of cuts, matrix of cuts, or vector of cuts that will be
#'    repeated.
#' @param type One of `"H"` and `"std"`, defaults to `"H"`.
#' @return The population ordinal reliabiltiy
#' @name ordinal_omega
#' @export

ordinal_omega = function(lambda, sigma, cuts, type = c("H", "std")) {

  checkmate::assertNumeric(lambda)
  k = length(lambda)
  checkmate::assertNumeric(sigma, len = k)
  cuts = massage_cuts(cuts, k)
  checkmate::assertList(cuts, len = k)
  type = match.arg(type)

  lambda = standardize_lambda(lambda, sigma)
  sigma = standardize_sigma(lambda, sigma)
  rho = tcrossprod(lambda, lambda) + diag(sigma^2)
  v = thurstone(lambda, sigma)

  xi = xi_theoretical(cuts, rho)

  if(type == "H") {
    c(crossprod(v, xi %*% v))
  } else {
    i = rep(1, k)
    c(crossprod(i, xi %*% v))^2/sum(xi)
  }

}

#' @rdname ordinal_omega
ordinal_std = function(lambda, sigma, cuts)
  ordinal_omega(lambda, sigma, cuts, type = "std")

#' @rdname ordinal_omega
ordinal_H = function(lambda, sigma, cuts)
  ordinal_omega(lambda, sigma, cuts, type = "H")

#' Empirical ordinal omega
#'
#' @param poly A an object of class `poly`.
#' @param y Optional array of observations.
#' @param xi_type One of `"theoretical"` and `"sample"`, defaults to
#'    `"theoretical"`. The option `"sample"` needs `y`
#' @param type One of `"H"` and `"std"`, defaults to `"H"`.
#' @param ... Passed to `psych::fa`.
#' @export
#' @return The empirical ordinal reliability.

ordinal_poly = function(poly, y, xi_type = c("theoretical", "sample"),
                        type = c("H", "std"), ...) {

  type = match.arg(type)
  xi_type = match.arg(xi_type)
  cuts = poly$tau
  rho = poly$rho
  k = nrow(rho)
  fa = psych::fa(poly$rho, ...)
  v = thurstone(fa$loadings, sqrt(fa$uniquenesses))

  if(xi_type  == "sample") {
    xi = xi_sample(y, cuts)
  } else if (xi_type  == "theoretical") {
    xi = xi_theoretical(cuts, rho)
  }

  if(type == "H") {
    c(crossprod(v, xi %*% v))
  } else {
    i = rep(1, k)
    c(crossprod(i, xi %*% v))^2/sum(xi)
  }

}

#' Calculating the Empirical Ordinal Alpha
#'
#' Use `ordinal_alpha` to calculate the empricial ordinal alpha from a `poly`
#'    object. The `poly` object must contain a named element `rho`, the
#'    polychoric correlation matrix, and a matrix `tau` with the associated
#'    cut points. See `[pscyh::polychoric]`.
#'
#' The population value of ordinal alpha equals the ordinal reliability when
#'    the underlying multivariate normal is parallel. The ordinal reliability
#'    is the sqaured correlation between the true latent variable and the
#'    best linear predictor of the observed Likert-type data.
#'
#' @param poly A an object of class `poly`.
#' @export
#' @return The empirical ordinal alpha.
#' @examples
#'  agreeableness = psychTools::bfi[c("A1", "A2", "A3", "A4", "A5")]
#'  agreeableness[, "A1"] = 7 - agreeableness[, "A1"] # Reverse-coded item.
#'  poly = psych::polychoric(agreeableness)
#'  ordinal_alpha(poly) # 0.6267724
#'  ordinal_poly(poly, type = "std") # 0.6394087

ordinal_alpha = function(poly) {

  cuts = poly$tau
  rho = poly$rho
  alpha_std(rho) * sum(xi_theoretical(cuts, rho)) / sum(rho)

}
