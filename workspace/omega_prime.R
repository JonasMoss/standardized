x_hat = function(y, cuts) {

  stopifnot(is.null(dim(y)))
  f = function(i) {
    -(dnorm(cuts[i + 1]) - dnorm(cuts[i]))/(pnorm(cuts[i + 1]) - pnorm(cuts[i]))
  }
  sapply(y, f)
}

x_hats = function(y, cuts) {
  apply(y, 2, x_hat, cuts)
}

omega_prime_hat = function(y, cuts) {
  xhats = x_hats(y, cuts)
  cov(xhats)
}

omega_prime_hat = function(y) {
  poly = psych::polychoric(y)
  rho = poly$rho
  k = ncol(y)
  n = nrow(y)
  cuts = cbind(rep(-Inf, k), poly$tau, rep(Inf, k))
  xhats = do.call(what = cbind,
            args = lapply(seq.int(k), function(i) x_hat(y[, i], cuts[i, ])))
  Xi = cov(xhats)

  lambda = psych::fa(rho)$loadings
  sigma = sqrt(psych::fa(rho)$uniquenesses)
  # tcrossprod(lambda, lambda) + diag(sigma^2)

  # This one find alpha prime.
  w0 = mean(lambda)/(k*mean(lambda)^2 + mean(sigma^2))
  sum(Xi) * w0^2

  w = lambda / (sigma^2 * (1 + sum(lambda * 1/sigma^2)))
  crossprod(rep(1, k), Xi %*% w)^2 / sum(Xi)
}

omega_prime(lambda, sigma, tau) {
  # Find w0
  lambda_star = lambda/sqrt(lambda^2 + sigma^2)
  sigma_star = sigma/sqrt(lambda^2 + sigma^2)
  w0 = mean(lambda)/(k*mean(lambda)^2 + mean(sigma^2))

}




x_hat(c(y[, 1]), cuts)
x_hats(y, cuts)

corr = cov(xhats)
colnames(corr) <- c("X1", "X2", "X3", "X4", "X5")
rownames(corr) <- colnames(corr)
model = " y =~ X1 + X2 + X3 + X4 + X5"

fit <- lavaan::sem(model = model, sample.cov = corr, sample.nobs = n)
coefs <- lavaan::lavInspect(fit, what = "x")

lambda <- abs(c(coefs$lambda * sqrt(as.numeric(coefs$psi))))
sigma <- sqrt(diag(lavaan::lavInspect(fit, what = "x")$theta))
