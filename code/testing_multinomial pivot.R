D <- matrix(rnorm(10000, 0, 1), ncol = 2)
D <- cbind(D, 0)

D_exp = exp(D)

D_norm = D_exp/rowSums(D_exp)

colMeans(D_norm)

apply(D_norm, 2, sd)
