simulation = function(nreps) {
  
  ##  The parameters used; just as in the paper!
  params = expand.grid(k = c(4, 20), n = c(20, 200), 
                       sigma = c(1, 0.5, 0.1),
                       error = c(error_t, error_beta, error_gamma))
  
  sims = t(apply(params, 1, function(param) {
    simulate_MSE(nreps = nreps, n = param$n, k = param$k, sigma = param$sigma,
                 error = param$error)
  }))
  
  ## This is a convoluted way of changing the sims result into the format used
  ## in the table in the paper. Please forgive my sloppiness!
  
  sims_1 = sims[, 1]
  dim(sims_1) = c(4, 9)
  sims_2 = sims[, 2]
  dim(sims_2) = c(4, 9)
  dim(sims) = c(8, 9)
  sims[c(1, 3, 5, 7), ] = sims_1
  sims[c(2, 4, 6, 8), ] = sims_2
  
  sims
  
}

set.seed(313)
simulation(1000)