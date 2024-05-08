
for(K in c(3,5,10)){
  
  D <- matrix(rnorm(10000*(K-1), 0, 3), ncol = K-1)
  D <- cbind(D, 0)
  
  D_exp = exp(D)
  
  D_norm = D_exp/rowSums(D_exp)
  
  cc = colMeans(D_norm)
  sds = apply(D_norm, 2, sd)
  
  print(sprintf('mean k: %s, mean K: %s', mean(cc[1:K-1]), cc[K]))
  print(sprintf('sd k: %s, sd K: %s', mean(sds[1:K-1]), sds[K]))
}



## Next simulation
# Simulate data WITHOUT pivot. So three independent N(0,1)s. Then fit a multinomial model on these with just the intercept. What are the results. Biased? Unbiased? I could set n to some number, like 50, and then run R = 1000 model fits. Does that make sense?

require(foreign)
require(nnet)
require(ggplot2)
require(reshape2)

n = 50
R = 1000

K = 5
res <- t(sapply(1:R, function(i){
  y <- sample(K, n, replace = T)
  mod_fit <- multinom(y ~ 1)
  if(mod_fit$convergence != 0){
    stop
  }
  coefs <- summary(mod_fit)$coefficients[,1]
  p_coefs <- c(1, exp(coefs))/(1 + sum(exp(coefs)))
  p_coefs
}))

colMeans(res)
apply(res, 2, sd)
# doesn't seem to be an issue here. What's wrong with the way I am thinking about this? Or is there a fancier model fit going on here?

# Venables, W. N. and Ripley, B. D. (2002) Modern Applied Statistics with S. Fourth edition. Springer. 

# That's where this comes from...


## Take 2: mlogit
library(mlogit)

D$y = factor(y, levels = as.character(1:10), ordered = T)

D2 = dfidx(D, shape = 'wide', choice = 'y')

tt <- summary(mlogit(y ~ 1, D2, reflevel = '10'))

res2 <- t(sapply(1:R, function(i){
  y <- sample(K, n, replace = T)
  D$y = factor(y, levels = as.character(1:K), ordered = T)
  D2 = dfidx(D, shape = 'wide', choice = 'y')
  mod_fit <- mlogit(y ~ 1, D2, reflevel = '1')

  coefs <- summary(mod_fit)$coefficients
  p_coefs <- c(1, exp(coefs))/(1 + sum(exp(coefs)))
  p_coefs
}))

# hmm no issue here. So the issue seems to be from simulating data rather than fitting. At least with how it's fit here. Maybe I could look into the notes for this
