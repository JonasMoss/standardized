source("functions.R")

plot_errors = function() {
  pdf("chunks/beta.pdf", height = 8, width = 8)

  par(las = 1)

  # sqrt(4*vbeta(1/10, 1/10))
  a = 1/10
  const = sqrt(vbeta(a, a))
  x = seq(-0.5/const, 0.5/const, by = 0.0001)
  plot(x = x,
       y = dbeta((x + 0.5/const)*const, a, a) * const,
       type = "l",
       cex.axis = 1.2,
       ylim = c(0, 24),
       lwd = 2,
       cex.lab = 1.2,
       xlab = "x",
       ylab = "Density")

  dev.off()
  pdf("chunks/gamma.pdf", height = 8, width = 8)

  a = 1/100
  const = sqrt(vgamma(a, a))
  x = seq(-1, const * 4, by = 0.001)/const
  plot(x = x,
       y = dgamma(x*const + 1, a, a) * const,
       type = "l",
       xlim = c(-1/const, 4),
       ylim = c(0, 10),
       cex.axis = 1.2,
       lwd = 2,
       cex.lab = 1.2,
       xlab = "x",
       ylab = "Density")
  dev.off()

}
plot_errors()


plot_ordinals = function() {
  pdf("chunks/ordinals.pdf", height = 8, width = 8)

  agreeableness = psychTools::bfi[c("A1", "A2", "A3", "A4", "A5")]
  agreeableness[, "A1"] = 7 - agreeableness[, "A1"] # Reverse-coded item.
  poly = psych::polychoric(agreeableness)
  fa = psych::fa(poly$rho,)
  lambda = c(fa$loadings)
  sigma = c(sqrt(fa$uniquenesses))

  outs = c(2:10, 15, 20, 25, 30, 40, 50)
  ordinal_Hs = sapply(outs, function(out) {
    ordinal_H(lambda, sigma, cuts = qnorm(seq(0, 1, length.out = out + 1)))
  })

  ordinal_omegas = sapply(outs, function(out) {
    ordinal_std(lambda, sigma, cuts = qnorm(seq(0, 1, length.out = out + 1)))
  })

  par(las = 1, mar = c(5.1, 4.1, 4.1, 4.1))
  plot(outs, ordinal_Hs, ylim = c(0.3, 0.8), pch = 19, log = "x",
       xlab = "Levels", ylab = "Reliability")
  points(outs, ordinal_omegas, pch = 1)
  abline(h = omega(lambda, sigma), lty = 2)
  abline(h = omega_h(lambda, sigma))
  labs = formatC(c(omega(lambda, sigma),
                   omega_h(lambda, sigma)), digits = 2)

  axis(side = 4, at = c(omega(lambda, sigma), omega_h(lambda, sigma)),
       labels = labs)

  dev.off()
}

plot_ordinals()
