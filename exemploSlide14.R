f=function(theta){
  a=0
  if(theta >2){a=1}
  
  return(a)
}

M=10000000

thetas = rcauchy(M, location = 0, scale = 1)
result = rep(NA, M)
for (i in 1:M) {
  result[i] = f(thetas[i])
}
mean(result)




# Modelo 2 (slide 14)
fq = function(u){
  ret = u^(-2)/(2*pi*(1+u^(-2)))
  return(ret)
}

M = 10000
us = runif(M, 0, 0.5)

result = rep(NA, M)
for (i in 1:M) {
  result[i] = fq(us[i])
}
mean(result)
