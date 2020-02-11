####################################################################################
#
#
# Bayesian Tobit regression (Sampling importance resampling)
#
#
# http://hedibert.org/wp-content/uploads/2013/12/limiteddependentvariablemodels.pdf
#
# Source:  1976 Panel Study of Income Dynamics,
# based on data for the previous year, 1975.  
# Of the 753 observations, the first 428 are for women with positive hours
# worked in 1975, while the remaining 325 observations are for women who
# did not work for pay in 1975.  A more complete discussion of the data is
# found in Mroz [1987], Appendix 1. 
#
# Thomas A. Mroz (1987) The Sensitivity of an Empirical Model of Married 
# Women's Hours of Work to Economic and Statistical Assumptions.  
# Econometrica, Vol. 55, No. 4 (Jul., 1987), pp. 765-799 
# Stable URL: http://www.jstor.org/stable/1911029
#
####################################################################################
#
# LFP  "A dummy variable = 1 if woman worked in 1975, else 0";
# WHRS "Wife's hours of work in 1975";
# KL6  "Number of children less than 6 years old in household";
# K618 "Number of children between ages 6 and 18 in household";
# WA   "Wife's age";
# WE   "Wife's educational attainment, in years";
# WW   "Wife's average hourly earnings, in 1975 dollars";
# RPWG "Wife's wage reported at the time of the 1976 interview
#       (not the same as the 1975 estimated wage).
#       To use the subsample with this wage, one needs to select 1975
#  workers with LFP=1, then select only those women with non-zero RPWG.
#  Only 325 women work in 1975 and have a non-zero RPWG in 1976.";
# HHRS "Husband's hours worked in 1975";
# HA   "Husband's age";
# HE   "Husband's educational attainment, in years";
# HW   "Husband's wage, in 1975 dollars";
# FAMINC "Family income, in 1975 dollars.\
#  This variable is used to construct the property income variable.";
# MTR  "This is the marginal tax rate facing the wife, and is taken from\
#  published federal tax tables (state and local income taxes are excluded).\
#  The taxable income on which this tax rate is calculated\
#  includes Social Security, if applicable to wife.";
# WMED "Wife's mother's educational attainment, in years";
# WFED "Wife's father's educational attainment, in years";
# UN   "Unemployment rate in county of residence, in percentage points.\
#  This taken from bracketed ranges.";
# CIT  "Dummy variable = 1 if live in large city (SMSA), else 0";
# AX   "Actual years of wife's previous labor market experience";
#
####################################################################################

# clear
rm(list = ls())

data = read.table("http://hedibert.org/wp-content/uploads/2020/01/mroz-data.txt",header=TRUE)
# attach(data)
# detach(data)
y = data$WHRS
x = data$AX

namex = "Years of previous experience"
namey = "Hours of work"

loglike = function(alpha, beta, sigma){
  sum(dnorm(y[y>0],alpha+beta*x[y>0],sigma,log=TRUE))+
    sum(pnorm(-(alpha+beta*x[y==0])/sigma,log=TRUE))
}

minusloglike = function(theta){
  alpha = theta[1]
  beta = theta[2]
  sigma = exp(theta[3])
  ret = -loglike(alpha, beta, sigma)
  
  return(ret)
  # -sum(dnorm(y[y>0],alpha+beta*x[y>0],sigma,log=TRUE))-
  #   sum(pnorm(-(alpha+beta*x[y==0])/sigma,log=TRUE))
}

theta = c(900, 30, log(800))
theta.mle = nlm(minusloglike, theta)$estimate
alpha.mle = theta.mle[1]
beta.mle = theta.mle[2]
sig.mle = exp(theta.mle[3])

par(mfrow=c(1,1))
plot(x,y,xlab=namex,ylab=namey)
abline(lm(y~x),col=2,lwd=2)
abline(lm(y[y>0]~x[y>0]),col=3,lwd=2)
points(x[y==0],y[y==0], pch=16)
points(x,alpha.mle+beta.mle*x,col=4, pch=16)
abline(alpha.mle, beta.mle,col=4,pch=16)
legend("topright",legend=c(
  paste("corr(x,y)=",round(cor(x,y),3),sep=""),
  paste("corr(x[y>0],y[y>0])=",round(cor(x[y>0],y[y>0]),3),sep=""),
  "Tobit likelihood"),col=2:4,lwd=2)

M = 10000
# Priori
#alphas = rnorm(M,-429.995052,10000)
#betas = rnorm(M,70.794740,1000)
#sigs = sqrt(1/rgamma(M,5/2,5*(1000^2)/2))

# Proposta
alphas = rnorm(M, alpha.mle, 1000)
betas = rnorm(M, beta.mle, 100)
sigs = sqrt(1/rgamma(M,5/2,5*(1000^2)/2))

# Inicializa vetor w
w = rep(0,M)


for (i in 1:M)
  w[i] = loglike(alphas[i],betas[i],sigs[i])+
  dnorm(alphas[i],-429.995052,10000,log=TRUE)+dnorm(betas[i],70.794740,1000,log=TRUE)-
  dnorm(alphas[i],-429.995052,1000, log=TRUE)-dnorm(betas[i],70.794740,100,log=TRUE)
w = exp(w-max(w))

# Amostragem do indice
ind = sample(1:M, size=M, replace=TRUE,prob=w)

# baseano no indice, busca-se alpha, beta e sigma
alphas1 = alphas[ind]
betas1 = betas[ind]
sigs1 = sigs[ind]

par(mfrow=c(1,2))
hist(alphas1,prob=TRUE)
lines(density(alphas),col=2)
hist(betas1,prob=TRUE)
lines(density(betas),col=2)


# AQUI DEVERIA SER ALPHAS1 E BETAS 1 NAO????
alpha.bayes = mean(alphas1)
beta.bayes = mean(betas1)


# x11()
par(mfrow=c(1,1))
plot(x,y,xlab=namex,ylab=namey)
abline(lm(y~x),col=2,lwd=2)
abline(lm(y[y>0]~x[y>0]),col=3,lwd=2)
points(x[y==0],y[y==0],pch=16)
legend("topright",legend=c(
  paste("corr(x,y)=",round(cor(x,y),3),sep=""),
  paste("corr(x[y>0],y[y>0])=",round(cor(x[y>0],y[y>0]),3),sep="")),col=2:3,lwd=2)
points(x,alpha.mle+beta.mle*x,col=4,pch=16)
points(x,alpha.bayes+beta.bayes*x,col=6,pch=16)
# points(x,mean(alphas1)+mean(betas1)*x,col=7,pch=16)


