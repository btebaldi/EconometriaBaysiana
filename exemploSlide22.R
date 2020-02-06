# Enveloping the standard normal density
# by a Cauchy density, or a uniform density

# Clear
rm(list=ls())

# Sampling importance resampling
# Objective: sample from N(0,1)
# Candidate: sample from U(-50,50)
# --------------------------------
# Bad proposal:U(-10, 10)
# Good proposal: Cauchy(0,1)
theta.grid = seq(from=-10, to=10, by=0.001)
N = length(theta.grid)
Au = max(dnorm(theta.grid)/(1/20))
1/Au
Ac = max(dnorm(theta.grid)/dcauchy(theta.grid))
1/Ac

# Limpeza de variaveis antigas
rm(list=ls())

