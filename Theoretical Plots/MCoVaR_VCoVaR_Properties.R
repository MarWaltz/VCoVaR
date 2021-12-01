#Note: The calculations in this script take only a couple of seconds/minutes on an Intel(R) 
#      Xeon(R) Gold 6136 CPU with 3.00 GHz on a 64-Bit Windows Server 2016 with 24 processors.

library("copula")

get_MultiCoVaR = function(copula, theta, alpha = 0.05, beta = 0.05){
  if(copula == "gumbel"){
    cop = gumbelCopula(param = theta, dim = 3)
  }
  if(copula == "clayton"){
    cop = claytonCopula(param = theta, dim = 3)
  }
  
  MinF = function(v){
    return(pCopula(c(v, alpha, alpha), copula = cop) - pCopula(c(1, alpha, alpha), copula = cop)*beta)
  }
  tmp = uniroot(f = MinF, lower = 0, upper = 1, tol = .Machine$double.eps*0.0000001, 
                maxiter = 10000000)$root
  MCoVaR = qnorm(tmp)
  return(MCoVaR)
}

get_VCoVaR = function(copula, theta, alpha = 0.05, beta = 0.05){
  if(copula == "gumbel"){
    cop = rotCopula(gumbelCopula(param = theta, dim = 3))
  }
  if(copula == "clayton"){
    cop = rotCopula(claytonCopula(param = theta, dim = 3))
  }
  
  MinF = function(v){
    return((v - pCopula(c(1, 1-alpha, 1-alpha), copula = cop) + 
              pCopula(c(1-v, 1-alpha, 1-alpha), copula = cop)) -
             (beta * (1 - pCopula(c(1, 1-alpha, 1-alpha), copula = cop))))
  }
  tmp = uniroot(f = MinF, lower = 0, upper = 1, tol = .Machine$double.eps*0.0000001, 
                maxiter = 10000000)$root
  VCoVaR = qnorm(tmp)
}

#Params
set.seed(1)
tau = seq(0.001, 0.999, length.out = 1000)

theta_Gu = iTau(copula = gumbelCopula(), tau = tau)
theta_Cl = iTau(copula = claytonCopula(), tau = tau)

### alpha = beta = 0.05 ###
#gumbel copula
M_gumbel05 = c()
for(i in seq_along(theta_Gu)){
  M_gumbel05[i] = get_MultiCoVaR(copula = "gumbel", theta = theta_Gu[i], 
                                 alpha = 0.05, beta = 0.05)
}
V_gumbel05 = c()
for(i in seq_along(theta_Gu)){
  V_gumbel05[i] = get_VCoVaR(copula = "gumbel", theta = theta_Gu[i], 
                             alpha = 0.05, beta = 0.05)
}

#clayton copula
M_clayton05 = c()
for(i in seq_along(theta_Cl)){
  M_clayton05[i] = get_MultiCoVaR(copula = "clayton", theta = theta_Cl[i], 
                                  alpha = 0.05, beta = 0.05)
}
V_clayton05 = c()
for(i in seq_along(theta_Cl)){
  V_clayton05[i] = get_VCoVaR(copula = "clayton", theta = theta_Cl[i], 
                              alpha = 0.05, beta = 0.05)
}


### alpha = beta = 0.01 ###
#gumbel copula
M_gumbel01 = c()
for(i in seq_along(theta_Gu)){
  M_gumbel01[i] = get_MultiCoVaR(copula = "gumbel", theta = theta_Gu[i], 
                                 alpha = 0.01, beta = 0.01)
}
V_gumbel01 = c()
for(i in seq_along(theta_Gu)){
  V_gumbel01[i] = get_VCoVaR(copula = "gumbel", theta = theta_Gu[i], 
                             alpha = 0.01, beta = 0.01)
}

#clayton copula
M_clayton01 = c()
for(i in seq_along(theta_Cl)){
  M_clayton01[i] = get_MultiCoVaR(copula = "clayton", theta = theta_Cl[i], 
                                  alpha = 0.01, beta = 0.01)
}
V_clayton01 = c()
for(i in seq_along(theta_Cl)){
  V_clayton01[i] = get_VCoVaR(copula = "clayton", theta = theta_Cl[i], 
                              alpha = 0.01, beta = 0.01)
}

#Remove MCoVaR values for Kendall's tau close to 1 (distorted from numerical issues)
M_gumbel01[970:1000] = NA
M_gumbel05[970:1000] = NA
M_clayton01[970:1000] = NA
M_clayton05[970:1000] = NA

#Comparison for 5% level
par(mfrow = c(1,2))
plot(tau, M_clayton05, type = "l", xlim = c(0, 1), ylab = "MCoVaR | VCoVaR", col = "red3", lwd = 2,
     xlab = expression(tau), lty = 2, ylim = c(-4, -1.5))
lines(tau, M_gumbel05, col = "blue3", lwd = 2, lty = 2)
lines(tau, V_clayton05, col = "red3", lwd = 2)
lines(tau, V_gumbel05, col = "blue3", lwd = 2)
legend("bottomright", col = rep(c("red3", "blue3"),2), bty = "n", density = rep(0,4), lwd = 2,
       legend = c("MCoVaR: clayton", "MCoVaR: gumbel", "VCoVaR: clayton", "VCoVaR: gumbel"), 
       border = c(NA,NA), lty = c(2,2,1,1))
mtext(text = expression(paste(alpha, " = ", beta, " = 0.05")))
abline(h = qnorm(0.05*0.05), lwd = 2, col = "darkgrey", lty = 2)
abline(h = qnorm(0.05), lwd = 2, col = "darkgrey", lty = 2)

#Comparison for 1% level
plot(tau, M_clayton01, type = "l", xlim = c(0, 1), ylab = "MCoVaR | VCoVaR", col = "red3", lwd = 2,
     xlab = expression(tau), lty = 2, ylim = c(-4, -1.5))
lines(tau, M_gumbel01, col = "blue3", lwd = 2, lty = 2)
lines(tau, V_clayton01, col = "red3", lwd = 2)
lines(tau, V_gumbel01, col = "blue3", lwd = 2)
legend("topright", col = rep(c("red3", "blue3"),2), bty = "n", density = rep(0,4), lwd = 2,
       legend = c("MCoVaR: clayton", "MCoVaR: gumbel", "VCoVaR: clayton", "VCoVaR: gumbel"), 
       border = c(NA,NA), lty = c(2,2,1,1))
mtext(text = expression(paste(alpha, " = ", beta, " = 0.01")))
abline(h = qnorm(0.01*0.01), lwd = 2, col = "darkgrey", lty = 2)
abline(h = qnorm(0.01), lwd = 2, col = "darkgrey", lty = 2)
