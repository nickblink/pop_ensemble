x <- (1:2000)/10
y1 <- dgamma(x = x, 1, 10)
y2 <- dgamma(x = x, 50, 0.5)
y3 <- dgamma(x = x, 1000, 10)
#y4 <- dgamma(x = x, .001, .001)

#pdf(file = 'my-gamma-plot.pdf')
plot(x, y1, type = 'l', ylim = c(0,0.2))
lines(x, y2, type = 'l', col = 2)
lines(x, y3, type = 'l', col = 3)
#lines(x, y4, type = 'l', col = 4)
legend("topright",
       legend = c('Gamma(1,10)', 'Gamma(50,0.5)', 'Gamma(1000,10)'),#, 'gamma4'), 
       col = c(1:4),
       lty = 1, cex = 0.8)


dgamma(x = 100, 1, 10)
dgamma(x = 100, 1000, 10)
