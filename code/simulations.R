Nreps = 10^5
set.seed(313)
sims = simulation(Nreps)

caption = "Simulations of $100 \\times \\textrm{MSE}_\\alpha/\\textrm{MSE}_{\\alpha_s}$ in the parallel model"

tab = prettify(sims * 100, 3, 0)

colnames(tab) = c("$t(3)$", "$t(3)$", "$t(3)$",
             "Beta", "Beta", "Beta",
             "Gamma", "Gamma", "Gamma")

rownames(tab) = c("$k = 5, n = 50$",
                  "$k = 20, n = 50$",
                  "$k = 5, n = 200$",
                  "$k = 20, n = 200$")

addtorow <- list()
addtorow$pos <- list(0, 0)
addtorow$command <- c(
  "Distribution & $t(5)$ & $t(5)$ & $t(5)$ & Beta & Beta & Beta & Gamma & Gamma & Gamma \\\\\n",
  "$\\sigma$ & $2$ & $1$ & $0.5$ & $2$ & $1$ & $0.5$ & $2$ & $1$ & $0.5$ \\\\\n"
)

tab = xtable::xtable(tab, caption = caption)
tab_str = print(tab,
      sanitize.colnames.function = identity,
      sanitize.names.function = identity,
      sanitize.rownames.function = identity,
      hline.after = NULL,
      caption.placement = "top",
      include.colnames = FALSE,
      add.to.row = addtorow,
      print.results = FALSE)

cat(tab_str, file = "chunks/simulations_table.tex")
cat(paste0("\\newcommand{\\geomean}{",round((prod(sims))^(1/length(sims)),2),"}"),
    file = "chunks/simulations_table.tex",
    append = TRUE,
    sep = "\n"
    )
