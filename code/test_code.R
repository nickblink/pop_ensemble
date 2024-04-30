

X <- rnorm(100)
beta0 = 2
beta1 = -1
y = beta0 + beta1*X + rnorm(100)

res <- lm(y ~ X)

save(res, file = 'results/test_04302024.RData')