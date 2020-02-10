# Exercise 6.1
# (C.) Jiazi Tang, 2009
rm(list=ls())
t=10^4
x[1]=rnorm(t)
r=0.9
for (i in 2:10^4){
  x[i]=r*x[i-1]+rnorm(1) 
  }
hist(x,freq=F,col="wheat2",main="")
curve(dnorm(x,sd=1/sqrt(1-r^2)),add=T,col="tomato")

rm(list = ls())
a = 2.7
b = 6.3
M = 200
x = rep(runif(1), M)

for (i in 2:M) {
  y = runif(1)
  rho = dbeta(y, a, b)/dbeta(x[i-1], a, b)
  rho = min(1,rho)
  x[i] = x[i-1] + (y-x[i-1])*(runif(1)<rho) 
}
x[1:5]
plot(x, type="l")


Nsim=5000
X=rep(runif(1),Nsim) # initialize the chain
for (i in 2:Nsim){
  Y=runif(1)
  rho=dbeta(Y,a,b)/dbeta(X[i-1],a,b)
  X[i]=X[i-1] + (Y-X[i-1])*(runif(1)<rho)
}
plot(X, type = "l")



rm(list = ls())
Nsim=10^4
X=c(rt(1,1)) # initialize the chain from the stationary
for (t in 2:Nsim){
  Y=rnorm(1) # candidate normal
  rho=dt(Y,1)*dnorm(X[t-1])/(dt(X[t-1],1)*dnorm(Y))
  X[t]=X[t-1] + (Y-X[t-1])*(runif(1)<rho)
}
plot(X[1:200], type = "l")

Nsim=10^4
X=c(rcauchy(1)) # initialize the chain from the stationary
for (t in 2:Nsim){
  Y    = rnorm(1) # candidate normal
  rho  = dcauchy(Y)*dnorm(X[t-1])/(dcauchy(X[t-1])*dnorm(Y))
  X[t] = X[t-1] + (Y-X[t-1])*(runif(1)<rho)
}
plot(X[1:200], type = "l")
