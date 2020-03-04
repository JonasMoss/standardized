# This is the code for the short simulations of section 3.

nreps = 10^5
k = 5

set.seed(313)
results1 = replicate(nreps, {
  lambda = rexp(k)
  sigma = rexp(k)
  c(omega_std = omega_std(lambda, sigma),
    omega_sigma = omega_sigma(lambda, sigma),
    omega = omega(lambda, sigma))
})

ranks1 = apply(-results1, 2, rank)

# How often is omega_sigma best?
mean(ranks1[2,] == 1) # [1] 0.98047
# How often is omega_std best?
mean(ranks1[1,] == 1) # [1] 0
# How often is omega best?
mean(ranks1[3,] == 1) # [1] 0.01953
# How often is omega larger than omega_std?
mean(results1[3, ] > results1[1, ]) # [1] 0.38878
# Mean reliability
rowMeans(results1)
# omega_std omega_sigma       omega
# 0.7644092   0.9091285   0.7003210

# Happens that omega_s < omega_std
sum(ranks1[2,] < ranks1[1,])


set.seed(313)
results2 = replicate(nreps, {
  lambda = runif(k)
  sigma = runif(k)
  c(omega_std = omega_std(lambda, sigma),
    omega_sigma = omega_sigma(lambda, sigma),
    omega = omega(lambda, sigma))
})

ranks2 = apply(-results2, 2, rank)

mean(ranks2[2,] == 1)
mean(ranks2[1,] == 1)
mean(ranks2[3,] == 1)
mean(results2[3, ] > results2[1, ])
rowMeans(results2)
sum(ranks2[2,] < ranks2[1,])
