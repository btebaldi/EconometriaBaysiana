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

# limpeza de variaveis
rm(list = ls())


# funcao de verossimilhanca
logden = function(lambda_i, y_i, x_i, beta, sigma)
{
  ret = -0.5*log(lambda_i) -0.5*(lambda_i+( (y_i - x_i %*% beta)/sigma  )^2 * 1/lambda_i)
  return(ret)
}



set.seed(12345)

n = 100
nregress = 3
lambda_true <- rexp(n,1/2)
sigma_true <- 2
error_true <- rnorm(n,0,sqrt(sigma_true*lambda_true))

X <- matrix(rnorm(3*n,1,2),100,nregress)
X[,1] <- 1

beta_true <- as.matrix(rep(runif(nregress),1))

Y <- X%*%beta_true + error_true  
plot(X%*%beta_true, Y)

ols <- lm(Y~X-1);
summary(ols)

hist(ols$residuals,breaks=15)

# priors
beta <- matrix(ols$coefficients,3,1)
sigma2 <- summary(ols)$sigma^2
sigma2_0 = 1

# MCMC set-up
M0    = 1000  # Final
M     = 10000 # Burn up
niter = M0+M

# Bayesian inference: 
# --------------------------------------
ncol = 1 + nregress + n
draws.ng = matrix(0, nrow = niter, ncol = ncol)
colnames(draws.ng) = c("sigma", paste("beta", 1:nregress), paste("Lambda", 1:n))
head(draws.ng)

# Bayesian inference: normal-gamma prior (using data augmentation)
# ----------------------------------------------------------------
# A SEREM DEFINIDOS
nu_0 = 2.5
V_0 = diag(sigma2, nregress)
beta_0 = beta

lambda = rep(1,n)

for (i in 1:(niter)){

  # Inicializo matrix Lambda^{-1}
  Lambda_1 = solve(diag(lambda))
  
  # full conditional of sigma2  
  d0=(nu_0 * sigma2_0)/2
  par1 = (nu_0 + n)/2
  par2 = d0 + ( t(Y-X%*%beta) %*% Lambda_1 %*% (Y-X%*%beta) )/2
  
  # Conditional distribution of sigma 
  sig2 = 1/rgamma(1,par1,par2)
  
  # full conditional of beta
  XtX = t(X) %*% Lambda_1 %*% X
  XtY = t(X) %*% Lambda_1 %*% Y
  
  # V_0 = diag(1/lambda)
  
  V_1   = solve(XtX/sig2 + solve(V_0))
  beta_1 = V_1 %*% (XtY + solve(V_0) %*% beta_0)
  
  beta =  beta_1 + t(chol(V_1)) %*% rnorm(nregress)
  
  # Seguindo terminologia de Koop. (Bayesian Econometrics Methods) pag 260
  
  # Draw of lambda
  for (j in 1:n){
    # psi1 = rexp(n, 0.5)
    y_j = Y[j, 1]
    x_j = X[j, ]
    
    # faco draw de no_0
    nu_michael_0 = rchisq(1,1)
    
    # calculo mu para linha j
    mu_j = abs(sqrt(sig2)/(y_j - x_j %*% beta)) 
    
    # Calculo x1 e x2
    x_1 = mu_j + (mu_j^2 * nu_michael_0)/2 - mu_j/2 * (4 * mu_j * nu_michael_0 + mu_j^2 * nu_michael_0^2)^0.5
    x_2 = mu_j^2 / x_1
    
    # decido quem sera escolhido
    p.treshold = mu_j/(mu_j + x_1)
    if (runif(1) < p.treshold){
      x_star = x_1
    } else{
      x_star = x_2
    }
    
    lambda_j = 1/x_star
    
    # draw from inverse Gausian
    if (lambda_j > 0){
      
      nume = logden(lambda_j,  y_j, x_j, beta, sqrt(sig2)) #% + dexp(lambda[j], rate = 0.5, log=TRUE)
      deno = logden(lambda[j], y_j, x_j, beta, sqrt(sig2)) #+ dexp(lambda_j , rate = 0.5, log=TRUE)
      
      log.rho = min(0, nume-deno)
      if (log(runif(1)) < log.rho){
        lambda[j] = lambda_j
      }
    }
  }
  
  # # storing draws
  draws.ng[i,] = c(sig2,beta, lambda)
}

draws2 = draws.ng[(M0+1):niter,]

summary(draws2)

# par(mfrow=c(2,2))
# ts.plot(draws1[,(p+2)],xlab="Iteration",ylab=expression(sigma),main="Without DA")
# acf(draws1[,(p+2)],main="")
# ts.plot(draws2[,(p+2)],xlab="Iteration",ylab=expression(sigma),main="With DA")
# acf(draws2[,(p+2)],main="")
# 
# par(mfrow=c(2,3))
# plot(density(draws1[,1]),main=expression(phi[0]),xlab="")
# lines(density(draws2[,1]),col=2)
# 
# plot(density(draws1[,2]),main=expression(phi[1]),xlab="")
# lines(density(draws2[,2]),col=2)
# 
# plot(density(draws1[,3]),main=expression(phi[2]),xlab="")
# lines(density(draws2[,3]),col=2)
# 
# plot(density(draws1[,4]),main=expression(phi[3]),xlab="")
# lines(density(draws2[,4]),col=2)
# 
# plot(density(draws1[,5]),main=expression(phi[4]),xlab="")
# lines(density(draws2[,5]),col=2)
# 
# plot(density(draws1[,p+2]),main=expression(sigma),xlab="")
# lines(density(draws2[,p+2]),col=2)