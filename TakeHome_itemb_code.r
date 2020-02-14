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

library(statmod)
library(nimble)

# limpeza de variaveis
rm(list = ls())


# funcao de verossimilhanca
func_q = function(lambda_i, y_i, x_i, beta, sigma)
{
  ret = -0.5*log(lambda_i) -0.5*(lambda_i+( (y_i - x_i %*% beta)/sigma  )^2 * 1/lambda_i)
  return(ret)
}

func_F = function(lambda_i, y_i, x_i, beta, sig2)
{
  erro = y_i - x_i%*%beta
  ret = dnorm(erro, mean = 0, sd = (sig2*lambda_i)^0.5, log = TRUE) - lambda_i/2 
  return(ret)
}

# Semente de randomizacao
set.seed(12345)

# tamanho da amostra
n = 200

# Numero de regressores
nregress = 4
beta_true = as.matrix(1:nregress, nregress,1)

# Regressors Matrix
X = matrix(rnorm(nregress*n,0,1), n, nregress)

# declare sigma
sigma_true = 1


# Simulacao
if(FALSE){
  lambda_true = rexp(n,1/2)
  error_true  = rnorm(n, 0, sqrt(sigma_true*lambda_true))
} else{
  error_true = nimble::rdexp(n, location = 0, scale = sig2^0.5)
}

Y <- X%*%beta_true + error_true  


# Ols de chute inicial
ols <- lm(Y~X-1);
summary(ols)

# histograma dos residuos
par(mfrow=c(1,2))
plot(X%*%beta_true, Y, xlab = "Y predicted", ylab = "Y")
hist(ols$residuals, breaks=15, main="Histogram of Ols resduals", xlab = "Residuals", ylab = "Frequency")

# PARAMETROS DE MCMC set-up
M0    = 1000  # Final
M     = 10000 # Burn up
niter = M0+M

# TABELA DE DRAWS: 
# --------------------------------------
ncol.draws = 1 + nregress
draws.ng = matrix(0, nrow = niter, ncol = ncol.draws)
colnames(draws.ng) = c("sigma", paste("beta", 1:nregress))

# priors of beta
beta_0 = matrix(0, nregress, 1)
V_0 = diag(25, nregress)

# priors of sigma
sigma2_0 = 1
nu_0 = 2.5

# initial Values
sigma2 = sigma2.ols
beta = beta.ols
lambda = rep(1,n)

for (i in 1:(niter)){

  # Inicializo matrix Lambda^{-1}
  Lambda_1 = solve(diag(lambda))
  
  # full conditional of sigma2  
  d0=(nu_0 * sigma2_0)/2
  par1 = (nu_0 + n)/2
  par2 = d0 + ( t(Y-X%*%beta) %*% Lambda_1 %*% (Y-X%*%beta) )/2
  
  # Conditional distribution of sigma 
  sig2 = 1/rgamma(1, par1, par2)
  
  # full conditional of beta
  XtX = t(X) %*% Lambda_1 %*% X
  XtY = t(X) %*% Lambda_1 %*% Y
  
  # V_0 = diag(1/lambda)
  
  V_1   = solve(XtX/sig2 + solve(V_0))
  beta_1 = V_1 %*% (XtY/sig2 + solve(V_0) %*% beta_0)
  
  beta =  beta_1 + t(chol(V_1)) %*% rnorm(nregress)
  
  # Seguindo terminologia de Koop. (Bayesian Econometrics Methods) pag 260
  
  # Draw of lambda
  for (j in 1:n){
    # psi1 = rexp(n, 0.5)
    y_j = Y[j, 1]
    x_j = X[j, ]
    
    # faco draw de no_0
    nu_michael_0 = rchisq(1,1)
    
    if(TRUE){
      #  x1 e x2
      x_1 = mu_j + (mu_j^2 * nu_michael_0)/2 - mu_j/2 * (4 * mu_j * nu_michael_0 + mu_j^2 * nu_michael_0^2)^0.5
      x_2 = mu_j^2 / x_1
      
      # decide between x_1 and x_2
      p.treshold = mu_j/(mu_j + x_1)
      if (runif(1) < p.treshold){
        x_star = x_1
      } else{
        x_star = x_2
      }
      
      # invert x_star
      lambda_j = 1/x_star
    }else
    {
      lambda_j = statmod::rinvgauss(1, mu_j, shape = 1) 
    }
    
    # draw from inverse Gausian
    if (lambda_j >0){
      
      # deno = func_F(lambda_j, y_j, x_j, beta, sig2) + func_q(lambda[j], y_j, x_j, beta, sqrt(sig2))
      # nume = func_q(lambda[j], y_j, x_j, beta, sqrt(sig2)) #+ dexp(lambda_j , rate = 0.5, log=TRUE)
      
      deno = func_F(lambda_j, y_j, x_j, beta, sig2) + func_q(lambda[j], y_j, x_j, beta, sqrt(sig2))
      nume = func_F(lambda[j], y_j, x_j, beta, sig2) + func_q(lambda_j, y_j, x_j, beta, sqrt(sig2))
      
      
      log.rho = min(0, nume-deno)
      if (log(runif(1)) < log.rho){
        lambda[j] = lambda_j
      }
    }
  }
  
  # storing draws
  draws.ng[i,] = c(sig2,beta, lambda)
}

# end_time <- 
Sys.time()
system.time()

draws2 = data.frame(draws.ng[(M0+1):niter,])
colnames(draws2)
summary(draws2$sigma)

x11()
par(mfrow=c(4,2))
hist(draws2$sigma, breaks = 50)
abline(v=sigma_true, col=2)
abline(v=sigma2, col="blue")
acf(draws2$sigma)

hist(draws2$beta.1, breaks = 50)
abline(v = beta_true[1], col=2)
abline(v = ols$coefficients[1], col="blue")
plot(draws2$beta.1,xlab="Iteration",ylab=expression(sigma),main="With DA", type = "l")

hist(draws2$beta.2, breaks = 50)
abline(v = beta_true[2], col=2)
abline(v = ols$coefficients[2], col="blue")
plot(draws2$beta.2,xlab="Iteration",ylab=expression(sigma),main="With DA", type = "l")

hist(draws2$beta.3, breaks = 50)
abline(v = beta_true[3], col=2)
abline(v = ols$coefficients[2], col="blue")
plot(draws2$beta.3,xlab="Iteration",ylab=expression(sigma),main="With DA", type = "l")
