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

dt_time = data.frame(Method = c("MH", "GIBBS"), Time = NA, ess=NA, ess_ps = NA)
for (k in 1:2) {
  if(k==1)
  { MH = TRUE }
  else
  { MH=FALSE }
  
  begin_time = Sys.time()
  for (i in 1:(niter)){
    
    # Inicialize Lambda^{-1}  matrix
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
    
    V_1   = solve(XtX/sig2 + solve(V_0))
    beta_1 = V_1 %*% (XtY/sig2 + solve(V_0) %*% beta_0)
    beta =  beta_1 + t(chol(V_1)) %*% rnorm(nregress)
    
    # Foolowing Koop. (Bayesian Econometrics Methods) pag 260
    
    # Draw of lambda
    for (j in 1:n){
      y_j = Y[j, 1]
      x_j = X[j, ]
      
      # draw de nu_0
      nu_michael_0 = rchisq(1,1)
      
      # mu for row j
      mu_j = abs(sqrt(sig2)/(y_j - x_j %*% beta)) 
      
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
      
      
      # now the Metropolis-Hasting
      if ((lambda_j >0) & (MH)){
        
        deno = func_F(lambda_j, y_j, x_j, beta, sig2) + func_q(lambda[j], y_j, x_j, beta, sqrt(sig2))
        nume = func_F(lambda[j], y_j, x_j, beta, sig2) + func_q(lambda_j, y_j, x_j, beta, sqrt(sig2))
        
        log.rho = min(0, nume-deno)
        if (log(runif(1)) < log.rho){
          lambda[j] = lambda_j
        }
      }
    }
    
    # storing draws
    if(MH)
    {draws.mc[i,] = c(sig2, beta)}
    else
    {draws.gibbs[i,] = c(sig2, beta)}
    
  }
  end_time = Sys.time()
  
  # Determine Execution Time
  dt_time[k, "Time"] = end_time - begin_time
  
  # Determine ESS
  dt_time[k, "ess"] = round(M/(1+2*sum(acf(draws.mc[,"sigma"],lag.max=1000,plot=FALSE)$acf[2:1001])))
}

# Determine the ess per second
dt_time$ess_ps = dt_time$ess / dt_time$Time

draws2 = data.frame(draws.mc[(M0+1):niter,])
colnames(draws2)
summary(draws2$sigma)

par(mfrow=c(1,2))
hist(draws2$sigma, breaks = 50, main="Histogram of Sigma", xlab = "Sigma", ylab = "Frequency", xlim = range(sigma_true, sigma2.ols, draws2$sigma))
abline(v=sigma_true, col="red")
abline(v=sigma2.ols, col="blue")
legend("bottom", c("OLS", "MCMC"), col = c("blue", "red"), lty=1)
acf(draws2$sigma)

par(mfrow=c(1,2))
hist(draws2$beta.1, breaks = 50, main="Histogram of beta_1", xlab = "beta_1", ylab = "Frequency")
abline(v = beta_true[1], col="red")
abline(v = beta.ols[1], col="blue")
legend("bottom", c("OLS", "MCMC"), col = c("blue", "red"), lty=1)
plot(draws2$beta.1, xlab="Iteration", ylab="beta_1",main="Interations", type = "l")

par(mfrow=c(1,2))
hist(draws2$beta.2, breaks = 50, main="Histogram of beta_2", xlab = "beta_2", ylab = "Frequency")
abline(v = beta_true[2], col="red")
abline(v = ols$coefficients[2], col="blue")
legend("bottom", c("OLS", "MCMC"), col = c("blue", "red"), lty=1)
plot(draws2$beta.2, xlab="Iteration", ylab="beta_2",main="Interations", type = "l")

par(mfrow=c(1,2))
hist(draws2$beta.3, breaks = 50, main="Histogram of beta_3", xlab = "beta_3", ylab = "Frequency")
abline(v = beta_true[3], col="red")
abline(v = ols$coefficients[3], col="blue")
legend("bottom", c("OLS", "MCMC"), col = c("blue", "red"), lty=1)
plot(draws2$beta.3, xlab="Iteration", ylab="beta_3",main="Interations", type = "l")

par(mfrow=c(1,2))
hist(draws2$beta.4, breaks = 50, main="Histogram of beta_4", xlab = "beta_4", ylab = "Frequency")
abline(v = beta_true[4], col="red")
abline(v = ols$coefficients[4], col="blue")
legend("bottom", c("OLS", "MCMC"), col = c("blue", "red"), lty=1)
plot(draws2$beta.3, xlab="Iteration", ylab="beta_4",main="Interations", type = "l")