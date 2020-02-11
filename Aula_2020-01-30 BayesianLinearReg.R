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
rm(list=ls())
den.ng = function(phi, lambda, gamma){
  abs(phi)^(gamma-0.5)*besselK(abs(phi)/lambda,gamma-1/2)/(sqrt(pi)*2^(gamma-0.5)*lambda^(gamma+0.5))
}

logden = function(psi, gamma, lambda, phi){
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
ols = lm(yy~X-1)
summary(ols)
phi.ols = ols$coef
sig.ols = summary(ols)$sigma
se.ols = sqrt(diag(solve(t(X)%*%X)))*sig.ols
qphi.ols = cbind(phi.ols+qnorm(0.025)*se.ols, phi.ols, phi.ols+qnorm(0.975)*se.ols)

par(mfrow=c(1,1))
plot(0:p,qphi.ols[,2],ylim=range(qphi.ols),pch=16,xlab="Lag",ylab="AR coefficient")
for (i in 1:(p+1))
  segments(i-1,qphi.ols[i,1],i-1,qphi.ols[i,3],lwd=2)
abline(h=0,lty=2)  

# Bayesian inference: hyperparameters of both priors
# --------------------------------------------------

# sigma
c0 = 2.5
d0 = 2.5
par1  = c0+(n-p)/2

# phi (normal-gamma prior)
lambda = 0.4
gamma  = 0.8

# phi (normal prior)
b0 = rep(0,p+1)
V0 = 2*gamma*lambda^2
iV0   = diag(1/V0,p+1)
b0iV0 = b0/V0

# Sufficient statistics 
XtX   = t(X)%*%X
Xty   = t(X)%*%yy

x11()
par(mfrow=c(1,1))
phi = seq(-3,3,length=1000)
plot(phi,den.ng(phi,lambda,gamma),type="l",lwd=2,xlab="",ylab="Density",main=expression(phi[j]))
lines(phi,dnorm(phi,b0,sqrt(V0)),col=2,lwd=2)
legend("topright",legend=c("Normal prior","Normal-Gamma prior"),col=2:1,lty=1,lwd=2)
points(ols$coef,rep(0,p+1),col=3,pch=16)


# MCMC set-up
M0    =  1000
M     = 10000
niter = M0+M

# Bayesian inference: normal prior
# --------------------------------
draws.n = matrix(0,niter,p+2)
phi     = rep(0,p+1)
for (i in 1:(niter)){
  # full conditional of sigma2
  par2 = d0+sum((yy-X%*%phi)^2)/2
  sig2 = 1/rgamma(1,par1,par2)
  
  # full conditional of phi
  V   = solve(XtX/sig2+iV0)
  m   = V%*%(Xty/sig2+b0iV0)
  phi = m + t(chol(V))%*%rnorm(p+1)
  
  # storing draws
  draws.n[i,] = c(phi,sig2)
}
draws.n = draws.n[(M0+1):niter,]
qphi.n = t(apply(draws.n[,1:(p+1)], 2, quantile, c(0.05,0.5,0.95)))

# Bayesian inference: normal-gamma prior
# --------------------------------------
sd.psi   = 0.01 
draws.ng = matrix(0,niter,p+2)
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
  draws.ng[i,] = c(phi,sig2)
}
draws.ng = draws.ng[(M0+1):niter,]
qphi.ng = t(apply(draws.ng[,1:(p+1)],2,quantile,c(0.025,0.5,0.975)))

par(mfrow=c(1,1))
limx = range(draws.n[,p+2],draws.ng[,p+2])
plot(density(draws.n[,p+2]),xlim=limx,main=expression(sigma),xlab="")
lines(density(draws.ng[,p+2]),col=2)
points(sig.ols,0,pch=16,col=3)

par(mfrow = c(1, 1))
plot(
  0:p,
  qphi.n[, 2],
  ylim = range(qphi.n, qphi.ng, qphi.ols),
  pch = 16,
  xlab = "Lag",
  ylab = "AR coefficient"
)
for (i in 1:(p + 1)) {
  segments(i - 1, qphi.n[i, 1], i - 1, qphi.n[i, 3], lwd = 2)
  segments(i - 1 + 0.25,
           qphi.ng[i, 1],
           i - 1 + 0.25,
           qphi.ng[i, 3],
           lwd = 2,
           col = 2)
  segments(i - 1 + 0.125,
           qphi.ols[i, 1],
           i - 1 + 0.125,
           qphi.ols[i, 3],
           lwd = 2,
           col = 3)
  points(i - 1 + 0.25, qphi.ng[i, 2], pch = 16, col = 2)
  points(i - 1 + 0.125, phi.ols[i], pch = 16, col = 3)
}
abline(h = 0, lty = 2)
legend(
  "topright",
  legend = c("Normal prior", "Normal-Gamma prior", "OLS"),
  col = 1:3,
  lty = 1,
  lwd = 2
)