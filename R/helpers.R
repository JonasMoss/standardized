#' Remove infinities from vector and append and prepend -Inf and Inf.
#'
#' @param x Numeric vector.
#' @return Fixed vector.

trim_vector = function(x) {
  checkmate::assertNumeric(x)
  x = x[x != Inf & x != -Inf]
  c(-Inf, x, Inf)
}

#' Massage cuts to the desired shape
#'
#' @param cuts A matrix, list, or vector of cuts
#' @param k Optional `k` saying how many times the vector of cuts should be
#'    repeated. Only matters when `cuts` is a vector.

massage_cuts = function(cuts, k)  {

  if (checkmate::testAtomicVector(cuts)) {
    if(missing(k)) k = 1
    cuts = t(matrix(rep(trim_vector(cuts), k), nrow = length(cuts)))
  }

  if(is.matrix(cuts)) {
    cuts = lapply(seq.int(nrow(cuts)), function(i) trim_vector(cuts[i, ]))
  }

  cuts

}

#' Calculate the theoretical Xi
#'
#' @param cuts A matrix, list, or vector of cuts
#' @param Sigma The polychoric correlation matrix.
#' @return Theoretical Xi matrix.

Xi_theoretical = function(cuts, Sigma) {

  items_k = nrow(Sigma)
  cuts = massage_cuts(cuts, items_k)
  probs = lapply(cuts, function(x) diff(stats::pnorm(x)))
  ms = sapply(probs, length)

  f = function(i, cut)
    -(stats::dnorm(cut[i + 1]) - stats::dnorm(cut[i])) / (stats::pnorm(cut[i + 1]) - stats::pnorm(cut[i]))

  etas = lapply(seq.int(probs), function(i) f(seq.int(ms[i]), cuts[[i]]))

  g = function(i, j, k, l) {
    lower = rep(-Inf, items_k)
    upper = rep(Inf, items_k)
    lower[c(i, j)] = c(cuts[[i]][k], cuts[[j]][l])
    upper[c(i, j)] = c(cuts[[i]][k + 1], cuts[[j]][l + 1])
    mvtnorm::pmvnorm(lower = lower, upper = upper, corr = Sigma)
  }

  FUN <- Vectorize(function(i, j) {

    if(i != j) {
      grid = as.matrix(expand.grid(seq.int(ms[i]), (seq.int(ms[j]))))
      sum(apply(grid, 1, function(y)
        g(i, j, y[1], y[2]) * etas[[i]][y[1]] * etas[[j]][y[2]]))
    } else {
      grid = seq.int(ms[i])
      sum(sapply(grid, function(y)
        g(i, i, y, y) * etas[[i]][y]^2))
    }

  })

  outer(X = seq.int(items_k), Y = seq.int(items_k), FUN = FUN)

}

#' Calculate sample XI
#'
#' @param y Array of observations.
#' @param cuts A matrix, list, or vector of cuts
#' @return Sample Xi matrix.

Xi_sample = function(y, cuts) {
  k = ncol(y)
  cuts = massage_cuts(cuts, k)
  args = lapply(seq.int(k), function(i) x_hat(y[, i], cuts[[i]]))
  xhats = do.call(what = cbind, args = args)
  stats::cov(xhats)
}

#' Transform Likert data to X_hats.
#'
#' @param y Vector of observations.
#' @param cuts Vector of cuts.
#' @return The X_hats associated with `y` and `cuts`.
x_hat = function(y, cuts) {

  stopifnot(is.null(dim(y)))

  f = function(i) {
    -(stats::dnorm(cuts[i + 1]) - stats::dnorm(cuts[i]))/
      (stats::pnorm(cuts[i + 1]) - stats::pnorm(cuts[i]))
  }

  sapply(y, f)

}

#' Standardize parameter vectors
#'
#' The function `standardize_lambda` standardizes `lambda` and
#'    `standardize_sigma` standardizes `sigma`.
#' @param lambda Vector of loadings.
#' @param sigma Vector of standard deviations.
#' @return Standardized vector.
#' @name standardize
standardize_lambda = function(lambda, sigma) lambda/sqrt(lambda^2 + sigma^2)

#' @rdname standardize
standardize_sigma = function(lambda, sigma) sigma/sqrt(lambda^2 + sigma^2)

