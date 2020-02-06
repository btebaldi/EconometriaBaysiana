# Monte Carlo approximaton 
# to E(theta) where theta ~ N(0,1)
# --------------------------------------------------
M = 100000
theta = rnorm(M)
mean(theta^2)

# Monte Carlo with Impotance Sampling 
# approximaton to E(theta) where theta ~ N(0,1)
# --------------------------------------------------

# q(theta) is Student t with 5 dof
theta1 = rt(M, 5)
mean(theta1^2)
mean(theta1^2*dnorm(theta1)/dt(theta1,5))

# q(theta) is Uniform(-100,100) 
theta2 = runif(M,-100, 100)
mean(theta2^2*dnorm(theta2)/(1/200))

# MC approximations to P(theta>2) where theta is Cauchy
# -----------------------------------------------------
M = 10000
theta = rt(M, 1)

I  = 1 - pt(2, 1)
I1 = cumsum(theta>2)/(1:M)
I2 = 0.5*(cumsum((theta<(-2)|(theta>2))))/(1:M)

plot(I1, ylim=range(I1,I2))
points(I2,col=2)
abline(h=I,col=4,lwd=3)


# Simple linear regression
# Using MCIS to compute posterior means
# -------------------------------------
set.seed(1245)
n = 100
x = runif(n)
y = rnorm(n,1+2*x,1)

plot(x,y)

lm(y~x)

post = function(alpha,beta){
  prod(dnorm(y,alpha+beta*x,1))*dnorm(alpha,0,10)*dnorm(beta,0,10)
}

M = 10000
as = rt(M,5)
bs = rt(M,5)
num1 = rep(0,M)
num2 = rep(0,M)
den = rep(0,M)
for (i in 1:M){
  num1[i] = as[i]*post(as[i],bs[i])/(dt(as[i],5)*dt(bs[i],5))
  num2[i] = bs[i]*post(as[i],bs[i])/(dt(as[i],5)*dt(bs[i],5))
  den[i] = post(as[i],bs[i])/(dt(as[i],5)*dt(bs[i],5))
}
alpha.p = mean(num1)/mean(den)
beta.p = mean(num2)/mean(den)


plot(x,y)
xnew = seq(0,1,length=100)
lines(xnew,alpha.p+beta.p*xnew,col=2)



# Sampling importance resampling
# Objective: sample from N(0,1)
# Candidate: sample from U(-50,50)
# --------------------------------

# Step 1: sample from candidate
M = 100000
theta = runif(M,-50,50)
hist(theta)

# Step 2: compute resampling weights
w = dnorm(theta)/(1/100)
hist(w)

# Step 3: resample
theta1 = sample(theta,replace=TRUE,size=1000,w)
hist(theta1)









