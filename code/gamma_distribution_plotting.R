x <- (1:200)/100
y1 <- dgamma(x = x, 1, 5)
y2 <- dgamma(x = x, 1, 1)
y3 <- dgamma(x = x, 2, 1)
y4 <- dgamma(x = x, .001, .001)

#pdf(file = 'my-gamma-plot.pdf')
plot(x, y1, type = 'l', ylim = c(0,1), main = 'pdf')
lines(x, y2, type = 'l', col = 2)
lines(x, y3, type = 'l', col = 3)
lines(x, y4, type = 'l', col = 4)
legend("topright",
       legend = c('Gamma(1,5)', 'Gamma(1,1)', 'Gamma(2,1)', 'Gamma(.001,.001)'), 
       col = c(1:4),
       lty = 1, cex = 0.8)

y1 <- pgamma(x, 1, 10)
y2 <- pgamma(x, 50, 0.5)
y3 <- pgamma(x, 1000, 10)
y4 <- pgamma(x, .001, .001)

plot(x, y1, type = 'l', ylim = c(0,1), main = 'cdf')
lines(x, y2, type = 'l', col = 2)
lines(x, y3, type = 'l', col = 3)
lines(x, y4, type = 'l', col = 4)
legend("bottomright",
       legend = c('Gamma(1,10)', 'Gamma(50,0.5)', 'Gamma(1000,10)', 'Gamma(.001,.001)'), 
       col = c(1:4),
       lty = 1, cex = 0.8)


dgamma(0.01, 1, 5)/dgamma(1, 1, 5)


dgamma(0.01, .002, .001)/dgamma(1, .002, .001)
