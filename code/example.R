#' Make pretty strings.
#' @param x Numeric matrix.
#' @param digits Digits to prettify with.

pretty = function(x, digits = 3, nsmall = 3) {
  y = sapply(x, function(x) format(signif(round(x, digits = digits), digits),
                                   nsmall = nsmall))
  dim(y) = dim(x)
  y
}


#' Make "x (y)" from x and y.
#' 
#' 
#' @param x Vector.
#' @param y Vector.
#' @value A vector of "mean (sd)".
#' @source https://stat.ethz.ch/pipermail/r-help/2011-July/284323.html

msd = function(means, sds) {
  stopifnot(dim(means) == dim(sds))
  y = paste(means, " (", sds, ")", sep = "") 
  dim(y) = dim(means)
  y
} 

models = list(' y  =~ A1 + A2 + A3 + A4 + A5 ',
              ' y  =~ C1 + C2 + C3 + C4 + C5 ',
              ' y  =~ E1 + E2 + E3 + E4 + E5 ',
              ' y  =~ N1 + N2 + N3 + N4 + N5 ',
              ' y  =~ O1 + O2 + O3 + O4 + O5 ')

reliabilities = sapply(models, function(model) {
  fit <- lavaan::cfa(model, data = psychTools::bfi)
  coefs <- lavaan::lavInspect(fit, what = "x")
  lambda <- c(coefs$lambda * sqrt(as.numeric(coefs$psi)))
  sigma <- sqrt(diag(lavaan::lavInspect(fit, what = "x")$theta))
  rel(lambda, sigma)[c(3, 4, 2, 1)]
})

tab = pretty(reliabilities)
tab = xtable::xtable(tab,
                     caption = "Comparison of reliability coefficients on personality data")

rownames(tab) = c("Congeneric reliability", 
                  "Coefficient \\it{H}", 
                  "Standardized reliability",
                  "Sigma reliability")

colnames(tab) = c("A", "C", "E", "N", "O")
xtable::align(tab) = c("l", rep("c", 5))
xtable::label(tab) = c("tab:reliabilites")
print(tab, sanitize.rownames.function = identity,
      sanitize.colnames.function = identity,
      hline.after = NULL,
      caption.placement = "top")
