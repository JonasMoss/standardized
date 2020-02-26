## This code can be used to generate examples where omega = alpha_std.

f = function(p, k = 3) {
  lambda = p[1:k]
  sigma = p[(k + 1):(2 * k)]
  Sigma = tcrossprod(lambda) + diag(sigma^2)
  (alpha_std(Sigma) - omega(lambda, sigma))^2
}

set.seed(313)
k = 3 # Number of tests
result = nlm(f, p = runif(2 * k), k = k)
lambda = result$estimate[1:k] # 0.45 0.19 0.73
sigma = result$estimate[(k + 1):(2 * k)] # 0.20 0.83 0.68
Sigma = tcrossprod(lambda) + diag(sigma^2)
alpha_std(Sigma) # 0.61
omega_std(lambda, sigma) # 0.61


