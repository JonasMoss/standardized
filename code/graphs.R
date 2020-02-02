source("functions.R")

plot_errors = function() {
  pdf("../plots/beta.pdf", height = 8, width = 8)
  
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
  pdf("../plots/gamma.pdf", height = 8, width = 8)
  
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
  
  # pdf("../plots/t.pdf", height = 8, width = 8)
  # x = seq(-3, 3, by = 0.001)
  # plot(x = x,
  #      y = dt(x, df = 3),
  #      type = "l",
  #      cex.axis = 1.2,
  #      lwd = 2,
  #      cex.lab = 1.2,
  #      xlab = "x",
  #      ylab = "Density")
  # dev.off()  
  
}
plot_errors()
