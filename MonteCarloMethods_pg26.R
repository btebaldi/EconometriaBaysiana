rm(list = ls())

#
# Two-component mixture of Gaussians (data augmentation)
#

# density of two-component mixture of Gaussians
dmixt = function(y, mu1, mu2){
  prod(p_mixt(y,mu1,mu2,1,1))
}

p_mixt = function(y, mu1, mu2, sd1, sd2){
  0.8*dnorm(y, mu1, sd1) + 0.2*dnorm(y, mu2, sd2)
}

# Mixture density
mu1=0
mu2=3
sd1=1
sd2=1

ys = seq(-5,15,length=1000)

par(mfrow=c(1,1))
plot(ys, p_mixt(ys, mu1, mu2, sd1, sd2),type="l",xlab="y",ylab="Density")


# Sampling from mixture of Gaussians
set.seed(12145)
n=20
mu = c(mu1,mu2)
sd = c(sd1,sd2)
z = rbinom(n, 1, 0.2)
y = rnorm(n,mu[z+1],sd[z+1])

par(mfrow=c(1,1))
hist(y,prob=TRUE,breaks=seq(-10,20,length=30),main="")
lines(ys,0.8*dnorm(ys,mu1,sd1)+0.2*dnorm(ys,mu2,sd2),col=2)

# Posterior of (mu1,mu2) (contour plot)
N = 100
mu1s = seq(-5,5,length=N)
mu2s = seq(-5,5,length=N)
post = matrix(0,N,N)
for (i in 1:N)
  for (j in 1:N)
    post[i,j] = dmixt(y,mu1s[i],mu2s[j])*dnorm(mu1s[i])*dnorm(mu2s[j])

contour(mu1s,mu2s,post,drawlabels=FALSE)



# Sampling importance resampling
M = 10000
mu1.d = rnorm(M)
mu2.d = rnorm(M)
w = rep(0,M)
for (i in 1:M)
  w[i] = dmixt(y,mu1.d[i],mu2.d[i])
ind = sample(1:M,size=M,replace=TRUE,prob=w)
mu1.d = mu1.d[ind]
mu2.d = mu2.d[ind]

par(mfrow=c(1,1))
plot(mu1.d,mu2.d)
contour(mu1s,mu2s,post,drawlabels=FALSE,add=TRUE,col=2,lwd=2)


# Running my first Gibbs Sampler
# Inicizlizo z (seletor de curva 0 ou 1?)
z = c(rep(0,n/2),rep(1,n/2))

# tamanho da amostra
M = 10000

# Inicializo a matrix de draws (m0, m1, z's)
draws = matrix(0,M,2+n)

# para cada linha
for (i in 1:M){
  
  # vejo quantos sao da curva 0
  n0 = sum(z==0)
  
  # media de y condicionado a cada curva
  ybar0 = sum(y[z==0])/(n0+1)
  ybar1 = sum(y[z==1])/(n-n0+1)
  
  # faz uma extracao da normal
  mu0 = rnorm(1,ybar0,sqrt(1/(n0+1)))
  mu1 = rnorm(1,ybar1,sqrt(1/(n-n0+1)))
  
  # calcula a probabilidade de cada amostra.
  # (veja que tenho de saber a 0.8 e 0.2)
  pz0 = dnorm(y,mu0,1)*0.8
  pz1 = dnorm(y,mu1,1)*0.2
  pz1 = pz1/(pz0+pz1)

  # Calculo novo z
  z   = rbinom(n, 1, pz1)
  
  # guardo a informação na matrix
  draws[i,] = c(mu0,mu1,z)
}

# determino limites de plotagem
limx = range(mu1.d,draws[,1])
limy = range(mu2.d,draws[,2])

par(mfrow=c(1,2))
plot(mu1.d,mu2.d,xlab=expression(mu[1]),ylab=expression(mu[2]),xlim=limx,ylim=limy)
title("SIR")
plot(draws[,1:2],xlab=expression(mu[1]),ylab=expression(mu[2]),xlim=limx,ylim=limy)
title("Gibbs sampler")









