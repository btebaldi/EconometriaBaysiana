########################################################################
#
#
# Our first Metropolis-Hastings algorithm
#
# Sampling from a bivariate 2-component mixture of Gaussians
#
########################################################################

rm(list = ls())

d2norm = function(theta){
  0.8*dnorm(theta[1])*dnorm(theta[2])+
    0.2*dnorm(theta[1], 1, 0.5)*dnorm(theta[2], 1, 0.5)
}

N = 100
thetas1 = seq(-5,5,length=N)
thetas2 = seq(-5,5,length=N)
den = matrix(0,N,N)
for (i in 1:N)
  for (j in 1:N)
    den[i,j] = d2norm(c(thetas1[i],thetas2[j]))

par(mfrow=c(1,1))
contour(den)

x11()
par(mfrow=c(4,3))
M = 1000
draws = rep(0,M)
for (sd in c(0.05,0.1,0.5,1)){
  contour(thetas1,thetas2,den)
  theta = c(4,-4)
  for (i in 1:M){
    
    # sorteio do theta* da distribuicao que tenho acesso
    theta.star = rnorm(2, theta, sd)
    
    # Numerador do rho
    nume=d2norm(theta.star)/prod(dnorm(theta.star,theta,sd))
    
    # Denominador do rho
    deno=d2norm(theta)/prod(dnorm(theta,theta.star,sd))
    
    # aka rho
    alpha = min(1,nume/deno)
    
    if (runif(1)<alpha){
      theta = theta.star
    }
    draws[i]=theta[1]
    points(theta[1],theta[2],col=3,pch=6)
  }
  # points(draws,theta[2],col=3,pch=16)
  contour(thetas1,thetas2,den,add=TRUE,col=2)
  title(paste("sd=",sd,sep=""))
  ts.plot(draws)
  acf(draws)
}




