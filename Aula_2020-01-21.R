rm(list = ls())

# Simulating some artificial data
set.seed(2325)
n   = 100
sig = 2
nu  = 4
y   = sig*rt(n,df=nu)
boxplot(y,horizontal=TRUE)


# SIR to draw from p(σ2|y,v)
# Below we implement SIR by sampling from a highly inefficient candidate density, U(0,100).

dt.hedi = function(sig2){
  prod(dt(y/sqrt(sig2),df=nu)/sqrt(sig2))
}

set.seed(4321)
M = 10000
sig2.t = runif(M,0,100)
w = rep(0,M)
for (i in 1:M)
  w[i] = dt.hedi(sig2.t[i])

sig2 = sample(sig2.t,size=M,replace=T,prob=w)

hist(sig2,xlab=expression(sigma^2),prob=TRUE)
abline(v=sig^2,col=2,lwd=3)




dt.hedi = function(sig2, nu){
  sum(dt(y/sqrt(sig2),df=nu,log=TRUE)-0.5*log(sig2))
}

M = 1000
sig2 = 1/rgamma(M,3.5,3.5)
nu.max = 100
nus = 1:nu.max
logpred = matrix(0,M,nu.max)
for (j in 1:M)
  for (i in 1:nu.max)
    logpred[j,i] = dt.hedi(sig2[j],nus[i])
A = max(logpred)
logpred1 = logpred-A
pred1 = exp(logpred1)
logpred1 = log(apply(pred1,2,sum)) + A -log(M)
plot(nus,logpred1,ylab="Log predictive",xlab=expression(nu))
abline(v=nu,col=2,lwd=3)