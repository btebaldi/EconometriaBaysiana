# Exemplo slide 9/40 Monte Carlo Method

# declaracao da funcao
f = function (theta)
{
  a=(cos(50*theta) + sin(20*theta))^2
  return(a)
}

# tamanho da amostra
M = seq(from=log(5), to=log(10000), by=0.05)
length(M)
M = unique(round(exp(M), digits=0))

convergence = rep(NA, length(M))
for (j in 1:length(M)) {

    # extracao da amostra
  thetas = runif(M[j], 0, 1)
  
  # inicializacao do vetor de retorno 
  # (calculo da funcao f para cada theta)
  result = rep(NA, M[j])
  
  for (i in 1:M[j])
  {
    result[i] = f(thetas[i])
  }
  
  convergence[j] = mean(result)
  
}

plot(M, convergence, type="l")
abline(a=0.965, b=0, col=2)

