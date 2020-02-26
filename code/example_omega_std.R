## Simulation from the second example.

## Preparations.
model <- ' y  =~ A1 + A2 + A3 + A4 + A5 '

fit <- lavaan::cfa(model, data = psychTools::bfi)
coefs <- lavaan::lavInspect(fit, what = "x")

lambda <- c(coefs$lambda * sqrt(as.numeric(coefs$psi)))
sigma <- sqrt(diag(lavaan::lavInspect(fit, what = "x")$theta))
omega_std_true <- omega_std(lambda, sigma)
alpha_std_true <- alpha_std(lambda %*% t(lambda) + diag(sigma^2))

## Simulation starts here.
set.seed(313)
n <- c(50, 200, 1000, 5000)
nreps = 1000

sims = sapply(n, function(n) {
  replicate(nreps, {
    sim <- simulate_tau(n, 5, lambda, sigma)
    colnames(sim) <- c("A1", "A2", "A3", "A4", "A5")
    fit <- lavaan::cfa(model, data = sim)
    coefs <- lavaan::lavInspect(fit, what = "x")
    lambda <- c(coefs$lambda * sqrt(as.numeric(coefs$psi)))
    sigma <- sqrt(diag(lavaan::lavInspect(fit, what = "x")$theta))
    c(omega_std = omega_std(lambda, sigma), alpha_std = alpha_std(cov(sim)))
  })
})

omegas = sims[2*(1:nreps) - 1, ]
alphas = sims[2*(1:nreps), ]

tab1 = rbind(n*colMeans((omegas - omega_std_true)^2, na.rm = TRUE),
             n*colMeans((alphas - omega_std_true)^2, na.rm = TRUE))
tab2 = rbind(colMeans((omegas - omega_std_true), na.rm = TRUE),
             colMeans((alphas - omega_std_true), na.rm = TRUE))

tab = rbind(formatC(tab1, digits = 2, drop0trailing = TRUE, format = "fg"),
            formatC(tab2[1, ], digits = 2, drop0trailing = TRUE, format = "g"),
            formatC(tab2[2, ], digits = 2, drop0trailing = TRUE, format = "fg"))
colnames(tab) = c("$n = 50$", "$n = 200$", "$n = 1000$", "$n = 5000$")
rownames(tab) = c("omega", "alpha", "omega", "alpha")
tab = print(xtable::xtable(tab), include.colnames = FALSE, hline.after = NULL,
      sanitize.rownames.function = identity,
      sanitize.colnames.function = identity,
      sanitize.text.function = identity,
      add.to.row = list(pos = list(0, 0, 2, 2),
      command = c("\\multicolumn{5}{c}{Mean squared error times $n$}\\tabularnewline \n",
                  "& $n = 50$ & $n = 200$ & $n = 1000$ & $n = 5000$ \\\\ \n ",
                  "\\multicolumn{5}{c}{Bias}\\tabularnewline \\\\ \n",
                  "& $n = 50$ & $n = 200\ & $n = 1000$ & $n = 5000$ \n ")))

tab = gsub("alpha.1", "$\\alpha_s$", tab, fixed = TRUE)
tab = gsub("omega.1", "$\\omega_s$", tab, fixed = TRUE)
tab = gsub("alpha", "$\\alpha_s$", tab, fixed = TRUE)
tab = gsub("omega", "$\\omega_s$", tab, fixed = TRUE)

cat(tab, file = "chunks/example_table_omega_alpha.tex")
