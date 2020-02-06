# ---------------------------------------
# Comparing two models for Bernoulli data
# ---------------------------------------

# Simulating some data
# --------------------
set.seed(12345)
n = 50
x = rep(0,n)
alpha = 0.0
beta = 1.0
for (t in 2:n){
  eta = alpha+beta*x[t-1]
  thetat = 1/(1+exp(-eta))
  x[t] = rbinom(1,1,thetat)
}
x11()
plot(x,xlab="Observation",ylab="")
title(paste("Bernoulli trials \n Fraction of 1s=",mean(x),sep=""))

# Likelihood function for model 1
# -------------------------------
M = 100
theta= seq(0, 1, length=M)
like.M1 = rep(0,M)
for (i in 1:M){
  like.M1[i] = prod(dbinom(x,1,theta[i]))
}
like.M1 = like.M1/max(like.M1)

par(mfrow=c(1,1))
plot(theta,like.M1,xlab=expression(theta),
     ylab="Likelihood",type="l")
title("Likelihood for model 1")

# Likelihood function for model 2
# -------------------------------
alphas = seq(-2,1,length=M)
betas = seq(0,3,length=M)
like.M2 = matrix(0,M,M)
like = rep(0, n-1)
for (i in 1:M){
  for (j in 1:M){
    for (t in 2:n){
      eta = alphas[i]+betas[j]*x[t-1]
      thetat = 1/(1+exp(-eta))
      like[t-1] = dbinom(x[t],1,thetat)
    }
    like.M2[i,j] = prod(like)
  }
}
like.M2 = like.M2/max(like.M2)

par(mfrow=c(1,2))
plot(theta,like.M1,xlab=expression(theta),
     ylab="Likelihood",type="l")
title("Likelihood for model 1")

contour(alphas,betas,like.M2,xlab=expression(alpha),
        ylab=expression(beta))
points(alpha,beta,pch=16,col=2)
title("Likelihood for model 2")


# MC-based inference (Model 1)
# ----------------------------
N = 50000
theta.draw = runif(N)
w = rep(0,N)
for (i in 1:N){
  w[i] = prod(dbinom(x, 1, theta.draw[i]))
}
ind = sample(1:N, size=N, replace=TRUE, prob=w)
theta.draw = theta.draw[ind]

par(mfrow=c(1,2))
plot(theta,like.M1,xlab=expression(theta),
     ylab="Likelihood",type="l")
title("Likelihood for model 1")
hist(theta.draw,prob=T,xlim=c(0,1),main="Posterior for model 1",xlab=expression(theta))


# MC-based posterior inference (Model 2)
# --------------------------------------
N = 50000
alpha.draw = runif(N,-4,3)
beta.draw  = runif(N,-2,5)
w = rep(0,N)
like = rep(0, n-1)
for (i in 1:N){
  for (t in 2:n){
    eta = alpha.draw[i]+beta.draw[i]*x[t-1]
    thetat = 1/(1+exp(-eta))
    like[t-1] = dbinom(x[t],1,thetat)
  }
  w[i] = prod(like)
}

ind = sample(1:N,size=N,replace=TRUE,prob=w)

alpha.draw = alpha.draw[ind]
beta.draw = beta.draw[ind]

# probabilidade dado que x_{t-1} = 0
eta0 = 1/(1+exp(-(alpha.draw)))
# probabilidade dado que x_{t-1} = 1
eta1 = 1/(1+exp(-(alpha.draw+beta.draw)))

par(mfrow=c(1,1))
plot(alpha.draw,beta.draw,xlab=expression(alpha),
     ylab=expression(beta))
contour(alphas,betas,like.M2,col=2,add=TRUE,lwd=3,drawlabels=T)
title("Likelihood & posterior for model 2")

# Posterior predictive for both models
# ------------------------------------
x11()
plot(density(eta0),xlim=c(0,1),ylim=c(0,7),xlab="Probability",ylab="Density",main="",lwd=2)
lines(density(eta1), col=2, lwd=2)
lines(density(theta.draw), col=4, lwd=2)
title("Pr(x[n+1] | x[1],...,x[n])")
legend("topleft",legend=c("Model 1","Model 2 (x[n]=0)","Model 2 (x[n]=1)"),col=c(4,1,2),lty=1,lwd=2)
points(mean(x),0,col=6,pch=16)
