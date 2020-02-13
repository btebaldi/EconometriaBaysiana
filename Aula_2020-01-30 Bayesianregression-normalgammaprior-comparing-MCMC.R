################################################################################
#
#  Bayesian linear regression
#  with Normal or Normal-Gamma prior for coefficients
#
################################################################################
#
# Analysing the Canadian Lynx data
# Moran (1953) fit an AR(2) model:
# y(t) = 1.05 + 1.41*y(t-1) - 0.77*y(t-2) + e(t)
#
# Moran (1953) The Statistical Analysis of the Canadian Lynx Cycle.
# Australian Journal of Zoology, 1, pp 163-173.
#
################################################################################
rm(list = ls())

den.ng = function(phi,lambda,gamma){
  abs(phi)^(gamma-0.5)*besselK(abs(phi)/lambda,gamma-1/2)/(sqrt(pi)*2^(gamma-0.5)*lambda^(gamma+0.5))
}

logden.phi = function(phi,m,v,gamma,lambda){
  dnorm(phi,m,sqrt(v),log=TRUE)+(gamma-0.5)*log(abs(phi))+
    log(besselK(abs(phi)/lambda,gamma-1/2))
}

logden = function(psi,gamma,lambda,phi){
  (gamma-1.5)*log(psi)-0.5*(psi/lambda^2+phi^2/psi)
}

data(lynx)
y = log10(lynx)
n = length(y)  
year = 1821:1934

par(mfrow=c(1,1))
plot(year,y,xlab="Year",ylab="Log number of trappings",main="",type="l")
title("Annual numbers of lynx trappings for 1821-1934 in Canada")

# Fitting an Gaussian AR(p) model
# -------------------------------
p  = 20
yy = y[(p+1):n]
X  = cbind(1,y[p:(n-1)])
for (k in 2:p)
  X = cbind(X,y[(p-k+1):(n-k)])

# ML estimation
phi.ols = lm(yy~X-1)$coef

# Bayesian inference: 
# Hyperparameter & sufficient statistics
# --------------------------------------
# sigma
c0 = 2.5
d0 = 2.5
par1  = c0+(n-p)/2

# phi
lambda = 0.4
gamma  = 0.8

# Sufficient statistics 
XtX   = t(X)%*%X
Xty   = t(X)%*%yy

# MCMC set-up
M0    = 1000
M     = 10000
niter = M0+M
sd.phi = 0.25

# Bayesian inference: normal-gamma prior
# --------------------------------------
sd.psi   = 0.01 
draws    = matrix(0,niter,p+2)
phi      = phi.ols
psi      = rep(1,p+1)
for (iter in 1:(niter)){
  # full conditional of sigma2
  par2 = d0+sum((yy-X%*%phi)^2)/2
  sig2 = 1/rgamma(1,par1,par2)
  
  # full conditional of phi
  V   = solve(XtX/sig2+diag(1/psi))
  m   = V%*%(Xty/sig2)
  phi = m + t(chol(V))%*%rnorm(p+1)
  
  for (i in 1:(p+1)){
    yyy = yy-X[,-i]%*%phi[-i]
    v = sig2/sum(X[,i]^2)
    m = sum(yyy*X[,i])/sum(X[,i]^2)
    phii = rnorm(1,phi[i],sd.phi)
    nume = logden.phi(phii,m,v,gamma,lambda)
    deno = logden.phi(phi[i],m,v,gamma,lambda)
    log.alpha = min(0,nume-deno)
    if (log(runif(1))<log.alpha){
      phi[i] = phii
    } 
  }
  
  # storing draws
  draws[iter,] = c(phi,sig2)
}
draws1 = draws[(M0+1):niter,]


# Bayesian inference: normal-gamma prior (using data augmentation)
# ----------------------------------------------------------------
sd.psi   = 0.01 
draws    = matrix(0,niter,p+2)
phi      = phi.ols
psi      = rep(1,p+1)
for (i in 1:(niter)){
  # full conditional of sigma2
  par2 = d0+sum((yy-X%*%phi)^2)/2
  sig2 = 1/rgamma(1,par1,par2)
  
  # full conditional of phi
  V   = solve(XtX/sig2+diag(1/psi))
  m   = V%*%(Xty/sig2)
  phi = m + t(chol(V))%*%rnorm(p+1)
  
  # full conditional of psi
  for (j in 1:(p+1)){
    psi1 = rlnorm(1,psi[j],sd.psi)
    if (psi1>0){
      nume = logden(psi1,gamma,lambda,phi[j])-dlnorm(psi1,psi[j],sd.psi,log=TRUE)
      deno = logden(psi[j],gamma,lambda,phi[j])-dlnorm(psi[j],psi1,sd.psi,log=TRUE)
      log.alpha = min(0,nume-deno)
      if (log(runif(1))<log.alpha){
        psi[j] = psi1
      }     
    } 
  }
  
  # storing draws
  draws[i,] = c(phi,sig2)
}
draws2 = draws[(M0+1):niter,]

par(mfrow=c(2,2))
ts.plot(draws1[,(p+2)],xlab="Iteration",ylab=expression(sigma),main="Without DA")
acf(draws1[,(p+2)],main="")
ts.plot(draws2[,(p+2)],xlab="Iteration",ylab=expression(sigma),main="With DA")
acf(draws2[,(p+2)],main="")

par(mfrow=c(2,3))
plot(density(draws1[,1]),main=expression(phi[0]),xlab="")
lines(density(draws2[,1]),col=2)

plot(density(draws1[,2]),main=expression(phi[1]),xlab="")
lines(density(draws2[,2]),col=2)

plot(density(draws1[,3]),main=expression(phi[2]),xlab="")
lines(density(draws2[,3]),col=2)

plot(density(draws1[,4]),main=expression(phi[3]),xlab="")
lines(density(draws2[,4]),col=2)

plot(density(draws1[,5]),main=expression(phi[4]),xlab="")
lines(density(draws2[,5]),col=2)

plot(density(draws1[,p+2]),main=expression(sigma),xlab="")
lines(density(draws2[,p+2]),col=2)