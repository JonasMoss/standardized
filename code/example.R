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
  c(omega(lambda, sigma),
    omega_h(lambda, sigma),
    omega_std(lambda, sigma),
    omega_sigma(lambda, sigma))
})

tab = prettify(reliabilities)
caption = "Comparison of reliability coefficients on personality data"
tab = xtable::xtable(tab, caption = caption)

rownames(tab) = c("Congeneric reliability",
                  "Coefficient \\it{H}",
                  "Standardized reliability",
                  "Sigma reliability")

colnames(tab) = c("A", "C", "E", "N", "O")
xtable::align(tab) = c("l", rep("c", 5))
xtable::label(tab) = c("tab:reliabilites")
tab_str = print(tab, sanitize.rownames.function = identity,
      sanitize.colnames.function = identity,
      hline.after = NULL,
      caption.placement = "top",
      print.results = TRUE)

description = "  \\vskip7.0pt
A, agreeableness; C, conscientiousness; E, extraversion; N, neuroticism; O, openness to experience \n"

tag = "\\end{tabular}\n"
strs = strsplit(tab_str, tag, fixed = TRUE)
cat(paste0(strs[[1]][1], tag, description, strs[[1]][2]),
    file = "chunks/example_table.tex")
