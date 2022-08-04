### Exploring the census data from 2010
### Will update this later to include Worldpop and FB (though I'll have to do a later year for FB at least)
library(ggplot2)

current_path <- getActiveDocumentContext()$path
setwd(dirname(current_path))

load('../pop_data/merged_denom_cov_data2206.RData')

# correlations
cor(adat[,3:5])

# looks like there is some bias
colMeans(adat[,3:5])

# plotting bias and such
adat$acs_bias = adat$acs - adat$census
adat$acs_prop_bas = adat$acs_bias/adat$census

adat$pep_bias = adat$pep - adat$census
adat$pep_prop_bas = adat$pep_bias/adat$census

adat$sd = apply(adat, 1, function(xx) sd(xx[3:5]))

# ggplot(data = adat) + 
#   geom_point(aes(x = census, y = acs_bias)) +
#   scale_x_continuous(trans='log2')


ggplot(data = adat) +
  geom_point(aes(x = census, y = abs(acs_bias))) +
  scale_x_continuous(trans='log2') +
  scale_y_continuous(trans='pseudo_log')

ggplot(data = adat) +
  geom_point(aes(x = census, y = abs(pep_bias))) +
  scale_x_continuous(trans='log2') +
  scale_y_continuous(trans='pseudo_log')

ggplot(data = adat) +  
  geom_point(aes(x = census, y = sd))  +
  scale_x_continuous(trans='log2') +
  scale_y_continuous(trans='log2')

