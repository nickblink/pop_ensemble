### Exploring the census data from 2010
### Will update this later to include Worldpop and FB (though I'll have to do a later year for FB at least)
library(ggplot2)

current_path <- rstudioapi::getActiveDocumentContext()$path
setwd(dirname(current_path))

load('../pop_data/merged_denom_cov_data2206.RData')

# correlations
cor(adat[,3:5])

# looks like there is some bias
colMeans(adat[,3:5])
colMeans(adat[,3:5])/mean(adat$census)

# plotting bias and such
adat$acs_bias = adat$acs - adat$census
adat$acs_prop_bas = adat$acs_bias/adat$census

adat$pep_bias = adat$pep - adat$census
adat$pep_prop_bas = adat$pep_bias/adat$census

adat$sd = apply(adat, 1, function(xx) sd(xx[3:5]))
adat$var = apply(adat, 1, function(xx) var(xx[3:5]))

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

ggplot(data = adat) +  
  geom_point(aes(x = census, y = var))  +
  geom_abline(intercept = 0, slope = 1, color = 'red') +
  scale_x_continuous(trans='log2') +
  scale_y_continuous(trans='log2') +
  ylab('Var(PEP, ACS, census)') +
  xlab('census count') +
  ggtitle('variance and census mean by county')

tt = adat %>%
  arrange(1-var)

acs_prop_bias <- sum(adat$acs)/sum(adat$census)
pep_prop_bias <- sum(adat$pep)/sum(adat$census)

adat_unbias <- adat %>%
  mutate(acs = acs/acs_prop_bias,
         pep = pep/pep_prop_bias)
adat_unbias$var = apply(adat_unbias, 1, function(xx) var(xx[3:5]))

ggplot(data = adat_unbias) +  
  geom_point(aes(x = census, y = var))  +
  geom_abline(intercept = 0, slope = 1, color = 'red') +
  scale_x_continuous(trans='log2') +
  scale_y_continuous(trans='log2') +
  ylab('Var(PEP, ACS, census)') +
  xlab('census count') +
  ggtitle('variance after bias adjusting')

tt = adat_unbias %>%
  arrange(1-var)
