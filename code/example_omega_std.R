## Simulation from the second example.

## Preparations.
model <- ' y  =~ A1 + A2 + A3 + A4 + A5 '

fit <- lavaan::cfa(model, data = psychTools::bfi)
coefs <- lavaan::lavInspect(fit, what = "x")

lambda <- abs(c(coefs$lambda * sqrt(as.numeric(coefs$psi))))
sigma <- sqrt(diag(lavaan::lavInspect(fit, what = "x")$theta))
omega_std_true <- omega_std(lambda, sigma)
psych::alpha(psychTools::bfi[, 1:5])
alpha_std_true <- alpha_std(lambda %*% t(lambda) + diag(sigma^2))

## Simulation starts here.
set.seed(313)
n <- c(50, 200, 1000, 5000)
nreps = 100

sims = sapply(n, function(n) {
  replicate(nreps, {
    sim <- simulate_tau(n, 5, lambda, sigma)
    colnames(sim) <- c("A1", "A2", "A3", "A4", "A5")
    fit <- lavaan::cfa(model, data = sim)
    coefs <- lavaan::lavInspect(fit, what = "x")
    lambda <- c(coefs$lambda * sqrt(as.numeric(coefs$psi)))
    sigma <- sqrt(diag(lavaan::lavInspect(fit, what = "x")$theta))
    c(omega_std = omega_std(lambda, sigma),
      alpha_std = alpha_std(cov(sim)))
  })
})

omegas = sims[2*(1:nreps) - 1, ]
alphas = sims[2*(1:nreps), ]

tab1 = rbind(n*colMeans((omegas - omega_std_true)^2, na.rm = TRUE),
             n*colMeans((alphas - omega_std_true)^2, na.rm = TRUE))
tab2 = rbind(colMeans((omegas - omega_std_true), na.rm = TRUE),
             colMeans((alphas - omega_std_true), na.rm = TRUE))

caption = "$n \\times \\textsc{MSE}$ for sample $\\omega_s$ and sample $\\alpha_s$"

tab = rbind(formatC(tab1, digits = 2, drop0trailing = TRUE, format = "fg"))
colnames(tab) = c("$n = 50$", "$n = 200$", "$n = 1000$", "$n = 5000$")
rownames(tab) = c("omegas", "alphas")
tab = print(xtable::xtable(tab, caption = caption, label = "tab:omega_std_alpha_std",),
      include.colnames = FALSE, hline.after = NULL,
      sanitize.rownames.function = identity,
      caption.placement = "top",
      sanitize.colnames.function = identity,
      sanitize.text.function = identity,
      add.to.row = list(pos = list(0),
      command = c("& $n = 50$ & $n = 200$ & $n = 1000$ & $n = 5000$ \\\\ \n ")),
      timestamp = paste0("Do not edit by hand. Last updated: ", date()))

tab = gsub("alphas", "$\\alpha_s$", tab, fixed = TRUE)
tab = gsub("omegas", "$\\omega_s$", tab, fixed = TRUE)

cat(tab, file = "chunks/example_table_omega_alpha.tex")
