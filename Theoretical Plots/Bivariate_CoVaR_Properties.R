#Note: The calculations in this script take only a couple of seconds/minutes on an Intel(R) 
#      Xeon(R) Gold 6136 CPU with 3.00 GHz on a 64-Bit Windows Server 2016 with 24 processors.

library("copula")

##################################### BIVARIATE COVAR #######################################
get_BivCoVaR = function(copula, theta, df = NULL, alpha = 0.05, beta = 0.05){
  if(copula == "normal"){
    MinF = function(v){
      return(pCopula(c(v, alpha), copula = normalCopula(param = theta, dispstr = "un")) - alpha*beta)
    }
  }
  if(copula == "t"){
    MinF = function(v){
      return(pCopula(c(v, alpha), copula = tCopula(param = theta, df = df, dispstr = "un")) - alpha*beta)
    }
  }
  if(copula == "gumbel"){
    MinF = function(v){
      return(pCopula(c(v, alpha), copula = gumbelCopula(param = theta)) - alpha*beta)
    }
  }
  if(copula == "clayton"){
    MinF = function(v){
      return(pCopula(c(v, alpha), copula = claytonCopula(param = theta)) - alpha*beta)
    }
  }
  tmp = uniroot(f = MinF, lower = 0, upper = 1, tol = .Machine$double.eps*0.0000001, maxiter = 10000)$root
  CoVaR = qnorm(tmp)
  
  return(CoVaR)
}

set.seed(1)
tau = seq(0,0.97,length.out = 1000)

theta_n = iTau(copula = normalCopula(), tau = tau)
theta_t = iTau(copula = tCopula(), tau = tau)
theta_Cl = iTau(copula = claytonCopula(), tau = tau)
theta_Gu = iTau(copula = gumbelCopula(), tau = tau)

### alpha = beta = 0.05 ###
#norm copula
Biv_norm05 = c()
for(i in seq_along(theta_n)){
  Biv_norm05[i] = get_BivCoVaR(copula = "normal", theta = theta_n[i], df = NULL,
                               alpha = 0.05, beta = 0.05)
}

#t copula (df = 4)
Biv_t05 = c()
for(i in seq_along(theta_t)){
  Biv_t05[i] = get_BivCoVaR(copula = "t", theta = theta_t[i], df = 4,
                            alpha = 0.05, beta = 0.05)
}

#t copula (df = 15)
Biv_t_205 = c()
for(i in seq_along(theta_t)){
  Biv_t_205[i] = get_BivCoVaR(copula = "t", theta = theta_t[i], df = 15,
                              alpha = 0.05, beta = 0.05)
}

#clayton copula
Biv_clayton05 = c()
for(i in seq_along(theta_Cl)){
  Biv_clayton05[i] = get_BivCoVaR(copula = "clayton", theta = theta_Cl[i], df = NULL,
                                  alpha = 0.05, beta = 0.05)
}

#gumbel copula
Biv_gumbel05 = c()
for(i in seq_along(theta_Gu)){
  Biv_gumbel05[i] = get_BivCoVaR(copula = "gumbel", theta = theta_Gu[i], df = NULL,
                                 alpha = 0.05, beta = 0.05)
}

### alpha = beta = 0.01 ###
#norm copula
Biv_norm01 = c()
for(i in seq_along(theta_n)){
  Biv_norm01[i] = get_BivCoVaR(copula = "normal", theta = theta_n[i], df = NULL,
                               alpha = 0.01, beta = 0.01)
}

#t copula (df = 4)
Biv_t01 = c()
for(i in seq_along(theta_t)){
  Biv_t01[i] = get_BivCoVaR(copula = "t", theta = theta_t[i], df = 4,
                            alpha = 0.01, beta = 0.01)
}

#t copula (df = 15)
Biv_t_201 = c()
for(i in seq_along(theta_t)){
  Biv_t_201[i] = get_BivCoVaR(copula = "t", theta = theta_t[i], df = 15,
                              alpha = 0.01, beta = 0.01)
}

#clayton copula
Biv_clayton01 = c()
for(i in seq_along(theta_Cl)){
  Biv_clayton01[i] = get_BivCoVaR(copula = "clayton", theta = theta_Cl[i], df = NULL,
                                  alpha = 0.01, beta = 0.01)
}

#gumbel copula
Biv_gumbel01 = c()
for(i in seq_along(theta_Gu)){
  Biv_gumbel01[i] = get_BivCoVaR(copula = "gumbel", theta = theta_Gu[i], df = NULL,
                                 alpha = 0.01, beta = 0.01)
}

#Comparison
pdf(file = "Sim_Bivariate_CoVaR.pdf", width = 12, height = 6)
par(mfrow = c(1,2))
plot(tau, Biv_norm05, type = "l", xlim = c(0, 1), ylab = "CoVaR", col = "black",
     lwd = 4, xlab = expression(tau), cex.axis = 1.1, cex.lab = 1.1)
lines(tau, Biv_t05, col = "red3", lwd = 4)
lines(tau, Biv_t_205, col = "blue3", lwd = 4)
lines(tau, Biv_gumbel05, col = "orange3", lwd = 4)
lines(tau, Biv_clayton05, col = "green3", lwd = 4)
abline(h = qnorm(0.05*0.05), lwd = 4, col = "darkgrey", lty = 2)
abline(h = qnorm(0.05), lwd = 4, col = "darkgrey", lty = 2)
legend("right", col = c("black", "red3", "blue3", "orange3", "green3"), bty = "n",
       legend = c("normal", "t (df = 4)", "t (df = 15)", "gumbel", "clayton"), 
       border = rep(NA,5), lty = rep(1,5), density = rep(0,5), lwd = 4)
mtext(text = expression(paste(alpha, " = ", beta, " = 0.05")), cex = 1.25)


plot(tau, Biv_norm01, type = "l", xlim = c(0, 1), ylab = "CoVaR", col = "black",
     lwd = 4, xlab = expression(tau), cex.axis = 1.1, cex.lab = 1.1)
lines(tau, Biv_t01, col = "red3", lwd = 4)
lines(tau, Biv_t_201, col = "blue3", lwd = 4)
lines(tau, Biv_gumbel01, col = "orange3", lwd = 4)
lines(tau, Biv_clayton01, col = "green3", lwd = 4)
abline(h = qnorm(0.01*0.01), lwd = 4, col = "darkgrey", lty = 2)
abline(h = qnorm(0.01), lwd = 4, col = "darkgrey", lty = 2)
legend("right", col = c("black", "red3", "blue3", "orange3", "green3"), bty = "n",
       legend = c("normal", "t (df = 4)", "t (df = 15)", "gumbel", "clayton"), 
       border = rep(NA,5), lty = rep(1,5), density = rep(0,5), lwd = 4)
mtext(text = expression(paste(alpha, " = ", beta, " = 0.01")), cex = 1.25)
dev.off()
