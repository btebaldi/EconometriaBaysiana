#########################################################
#
# SIR vs Gibbs sampler in the simple iid Gaussian case:
#
#           y1,...,yn iid N(mu,sig2)
#
# with prior:
#
#           p(mu,sig2) propto 1/sig2
#
#########################################################

# Clear
rm(list = ls())

# Simulating some Gaussian data
set.seed(54321)
n = 100
mu = 0
sig = 1
y = rnorm(n,mu,sig)
mu.mle = mean(y)
sig.mle = sqrt(var(y))


# SIR
sir.time = system.time({
  sumy1 = sum(y)
  sumy2 = sum(y^2)
  set.seed(12345)
  M0   = 600000
  M    = 10000
  mus  = runif(M0,-1,1)
  sigs = runif(M0,0,2)
  w = rep(0,M0)
  for (i in 1:M0)
    #  w[i] = sum(dnorm(y,mus[i],sigs[i],log=TRUE))-2*log(sigs[i])
    w[i] = -(n+2)*log(sigs[i])-0.5*(sumy2-2*mus[i]*sumy1+n*mus[i]^2)/(sigs[i]^2)
  w     = exp(w-max(w))
  ind   = sample(1:M0,size=M,prob=w,replace=TRUE)
  mus1  = mus[ind]
  sigs1 = sigs[ind]
})
sir.time = as.numeric(sir.time[3])
essps.Emu.sir = M/sir.time

par(mfrow=c(1,1))
plot(mus1,sigs1,col=2,pch=16,xlab=expression(mu),ylab=expression(sigma))

# Gibbs sampler
gibbs.time = system.time({
  n1   = (n-1)/2
  mu   = 0
  M0   = 10000
  M    = 10000
  draws = matrix(0,M0+M,2)
  for (i in 1:(M0+M)){
    sig2 = 1/rgamma(1,n1,sum((y-mu)^2)/2)
    mu   = rnorm(1,mu.mle,sqrt(sig2/n))
    draws[i,] = c(mu,sqrt(sig2))
  }
  draws = draws[(M0+1):(M0+M),]
})
gibbs.time = as.numeric(gibbs.time[3])
ess.Emu.gibbs = round(M/(1+2*sum(acf(draws[,1],lag.max=1000,plot=FALSE)$acf[2:1001])))
essps.Emu.gibbs = ess.Emu.gibbs/gibbs.time

table = cbind(c(sir.time,M,essps.Emu.sir),
              c(gibbs.time,ess.Emu.gibbs,essps.Emu.gibbs))

rownames(table) = c("Time (sec)","ESS","ESS/sec")
colnames(table) = c("SIR","Gibbs")
table

par(mfrow=c(1,2))
plot(mus1,sigs1,main="SIR",xlim=range(mus1,draws[,1]),ylim=range(sigs1,draws[,2]),
     xlab=expression(mu),ylab=expression(sigma))
plot(draws,main="Gibbs sampler",xlab=expression(mu),ylab=expression(sigma),
     xlim=range(mus1,draws[,1]),ylim=range(sigs1,draws[,2]))

x11()
par(mfrow=c(1,2))
plot(density(mus1),xlab="",ylab="Posterior density",main=expression(mu))
lines(density(draws[,1]),col=2)
legend("topleft",legend=c("SIR","Gibbs"),col=1:2,lty=1)
plot(density(sigs1),xlab="",ylab="Posterior density",main=expression(sigma))
lines(density(draws[,2]),col=2)

colMeans(tail(draws, 1000))