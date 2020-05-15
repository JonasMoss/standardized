##
## Only run the simulation if you have to.
##

if(!file.exists("chunks/simulations.Rds")) {
  Nreps = 100
  set.seed(313)
  simulations = simulation(Nreps)
  save(simulations, file = "chunks/simulations.Rds")
} else {
  load(file = "chunks/simulations.Rds")
}

##
## Make the table!
##

caption = "Simulations of $100 \\times \\textrm{MSE}_Z(\\alpha)/\\textrm{MSE}_Z(\\alpha_s)$ in the parallel model"

tab = prettify(simulations * 100, 3, 0)

colnames(tab) = c("$t(3)$", "$t(3)$", "$t(3)$",
             "Beta", "Beta", "Beta",
             "Gamma", "Gamma", "Gamma")

rownames(tab) = c("$k = 5, n = 50$",
                  "$k = 20, n = 50$",
                  "$k = 5, n = 200$",
                  "$k = 20, n = 200$")

addtorow <- list()
addtorow$pos <- list(0, 0, 0, 0, 4)
addtorow$command <- c(
  " \\toprule\n",
  " & $t(5)$ & $t(5)$ & $t(5)$ & Beta & Beta & Beta & Gamma & Gamma & Gamma \\\\\n",
  " \\cmidrule(lr){2-4}\\cmidrule(lr){5-7}\\cmidrule(lr){8-10} \n",
  " & $\\sigma = 2$ & $\\sigma = 1$ & $\\sigma = .5$ & $\\sigma = 2$ & $\\sigma = 1$ & $\\sigma = .5$ & $\\sigma = 2$ & $\\sigma = 1$ & $\\sigma = .5$ \\\\\n",
  " \\bottomrule\n"
)

tab = xtable::xtable(tab, caption = caption)

xtable::label(tab) = c("tab:simulation")

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
# cat(paste0("\\renewcommand{\\geomean}{$",
#            round((prod(simulations))^(1/length(simulations)),2),
#            "$}"),
#     file = "chunks/simulations_table.tex",
#     append = TRUE,
#     sep = "\n"
#     )
