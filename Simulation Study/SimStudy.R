library("copula")

eval = function(retX, retY, VaRX, CoVaR){
  exc = which(retX <= VaRX)
  
  red_CoVaR = rep(NA, length(CoVaR))
  red_CoVaR[exc] = CoVaR[exc]
  CoVaR_exc = which(retY <= red_CoVaR)
  CoVaR_rate = length(CoVaR_exc)/length(exc)
  
  return(CoVaR_rate)
}

evalMulti = function(retX1, retX2, retY, VaRX1, VaRX2, MCoVaR){
  exc1 = which(retX1 <= VaRX1)
  exc2 = which(retX2 <= VaRX2)
  exc = sort(intersect(exc1, exc2))
  
  red_MCoVaR = rep(NA, length(MCoVaR))
  red_MCoVaR[exc] = MCoVaR[exc]
  MCoVaR_exc = which(retY <= red_MCoVaR)
  MCoVaR_rate = length(MCoVaR_exc)/length(exc)
  
  return(MCoVaR_rate)
}

evalVulnerability = function(retX1, retX2, retY, VaRX1, VaRX2, VCoVaR){
  exc1 = which(retX1 <= VaRX1)
  exc2 = which(retX2 <= VaRX2)
  exc = sort(unique(c(exc1, exc2)))
  
  red_VCoVaR = rep(NA, length(VCoVaR))
  red_VCoVaR[exc] = VCoVaR[exc]
  VCoVaR_exc = which(retY <= red_VCoVaR)
  VCoVaR_rate = length(VCoVaR_exc)/length(exc)
  
  return(VCoVaR_rate)
}

OneStudyBiv = function(copula, tau, alpha, beta, sample.size = 10000, seed){
  if(copula == "clayton"){
    param = iTau(copula = claytonCopula(), tau = tau)  
  }
  if(copula == "gumbel"){
    param = iTau(copula = gumbelCopula(), tau = tau)
  }
  
  res = c()
  set.seed(seed)
  
  for(i in 1:100){
    if(copula == "clayton"){
      repeat{
        #Sample data
        sample = rCopula(n = sample.size, copula = claytonCopula(param = param))
        colnames(sample) = c("y", "x")
        
        #Estimate copula
        copEST = try(fitCopula(copula = claytonCopula(), data = sample, 
                               estimate.variance = F, method = "ml"), silent = T)
        if(class(copEST) != "try-error"){break}
      }
      cop = claytonCopula(param = copEST@estimate)
    }
    
    if(copula == "gumbel"){
      repeat{
        #Sample data
        sample = rCopula(n = sample.size, copula = gumbelCopula(param = param))
        colnames(sample) = c("y", "x")
        
        #Estimate copula
        copEST = try(fitCopula(copula = gumbelCopula(), data = sample, 
                               estimate.variance = F, method = "ml"), silent = T)
        if(class(copEST) != "try-error"){break}
      }
      cop = gumbelCopula(param = copEST@estimate)
    }
    
    #Calculate CoVaR and VaR
    MinF = function(v){
      return(pCopula(c(v, alpha), copula = cop) - alpha*beta)
    }
    CoVaR = uniroot(f = MinF, lower = 0, upper = 1, tol = .Machine$double.eps*0.0000001, 
                    maxiter = 10000000)$root
    VaRX = quantile(sample[,"x"], probs = alpha)
    
    #Evaluate
    tmp = eval(retX = sample[,"x"], retY = sample[,"y"], VaRX = VaRX, 
               CoVaR = rep(CoVaR, sample.size))
    res = c(res, tmp)
    cat(i)
  }
  
  return(mean(res))
}

OneStudyMulti = function(copula, tau, alpha, beta, sample.size = 10000, seed){
  if(copula == "clayton"){
    param = iTau(copula = claytonCopula(dim = 3), tau = tau)  
  }
  if(copula == "gumbel"){
    param = iTau(copula = gumbelCopula(dim = 3), tau = tau)
  }
  
  res = c()
  set.seed(seed)
  
  for(i in 1:100){
    if(copula == "clayton"){
      repeat{
        #Sample data
        sample = rCopula(n = sample.size, copula = claytonCopula(param = param, dim = 3))
        colnames(sample) = c("y", "x1", "x2")
        
        #Estimate copula
        copEST = try(fitCopula(copula = claytonCopula(dim = 3), data = sample, 
                               estimate.variance = F, method = "ml"), silent = T)
        if(class(copEST) != "try-error"){break}
      }
      cop = claytonCopula(param = copEST@estimate, dim = 3)
    }
    
    if(copula == "gumbel"){
      repeat{
        #Sample data
        sample = rCopula(n = sample.size, copula = gumbelCopula(param = param, dim = 3))
        colnames(sample) = c("y", "x1", "x2")
        
        #Estimate copula
        copEST = try(fitCopula(copula = gumbelCopula(dim = 3), data = sample, 
                               estimate.variance = F, method = "ml"), silent = T)
        if(class(copEST) != "try-error"){break}
      }
      cop = gumbelCopula(param = copEST@estimate, dim = 3)
    }
    
    #Calculate MCoVaR and VaR
    MinF = function(v){
      return(pCopula(c(v, alpha, alpha), copula = cop) - 
               pCopula(c(1, alpha, alpha), copula = cop)*beta)
    }
    MCoVaR = uniroot(f = MinF, lower = 0, upper = 1, tol = .Machine$double.eps*0.0000001, 
                     maxiter = 10000000)$root
    VaRX1 = quantile(sample[,"x1"], probs = alpha)
    VaRX2 = quantile(sample[,"x2"], probs = alpha)
    
    #Evaluate
    tmp = evalMulti(retX1 = sample[,"x1"], retX2 = sample[,"x2"], retY = sample[,"y"],
                    VaRX1 = VaRX1, VaRX2 = VaRX2, MCoVaR = rep(MCoVaR, sample.size))
    res = c(res, tmp)
    cat(i)
  }
  return(mean(res))
}

OneStudyVulnerability = function(copula, tau, alpha, beta, sample.size = 10000, seed){
  if(copula == "clayton"){
    param = iTau(copula = claytonCopula(dim = 3), tau = tau)  
  }
  if(copula == "gumbel"){
    param = iTau(copula = gumbelCopula(dim = 3), tau = tau)
  }
  
  res = c()
  set.seed(seed)
  
  for(i in 1:100){
    if(copula == "clayton"){
      repeat{
        #Sample data
        sample = rCopula(n = sample.size, copula = claytonCopula(param = param, dim = 3))
        colnames(sample) = c("y", "x1", "x2")
        
        #Estimate copula
        copEST = try(fitCopula(copula = claytonCopula(dim = 3), data = sample, 
                               estimate.variance = F, method = "ml"), silent = T)
        if(class(copEST) != "try-error"){break}
      }
      cop = rotCopula(claytonCopula(param = copEST@estimate, dim = 3))
    }
    
    if(copula == "gumbel"){
      repeat{
        #Sample data
        sample = rCopula(n = sample.size, copula = gumbelCopula(param = param, dim = 3))
        colnames(sample) = c("y", "x1", "x2")
        
        #Estimate copula
        copEST = try(fitCopula(copula = gumbelCopula(dim = 3), data = sample, 
                               estimate.variance = F, method = "ml"), silent = T)
        if(class(copEST) != "try-error"){break}
      }
      cop = rotCopula(gumbelCopula(param = copEST@estimate, dim = 3))
    }
    
    #Calculate VCoVaR and VaR
    MinF = function(v){
      return((v - pCopula(c(1, 1-alpha, 1-alpha), copula = cop) + 
                pCopula(c(1-v, 1-alpha, 1-alpha), copula = cop)) -
               (beta * (1 - pCopula(c(1, 1-alpha, 1-alpha), copula = cop))))
    }
    VCoVaR = uniroot(f = MinF, lower = 0, upper = 1, tol = .Machine$double.eps*0.0000001, 
                     maxiter = 10000000)$root
    VaRX1 = quantile(sample[,"x1"], probs = alpha)
    VaRX2 = quantile(sample[,"x2"], probs = alpha)
    
    #Evaluate
    tmp = evalVulnerability(retX1 = sample[,"x1"], retX2 = sample[,"x2"], retY = sample[,"y"],
                            VaRX1 = VaRX1, VaRX2 = VaRX2, VCoVaR = rep(VCoVaR, sample.size))
    res = c(res, tmp)
    cat(i)
  }
  return(mean(res))
}


#Bivariate CoVaR for alpha = beta = 0.05
BivTable_05 = matrix(NA, nrow = 2, ncol = 3, dimnames = list(c("clayton", "gumbel"),
                                                             c("0.25", "0.5", "0.75")))
for(i in seq_len(nrow(BivTable_05))){
  for(j in seq_len(ncol(BivTable_05))){
    copula = rownames(BivTable_05)[i]
    tau = as.numeric(colnames(BivTable_05)[j])
    BivTable_05[i,j] = OneStudyBiv(copula = copula, tau = tau, seed = 1 + i*j,
                                   alpha = 0.05, beta = 0.05)
  }
}

#Multi CoVaR for alpha = beta = 0.05
MultiTable_05 = matrix(NA, nrow = 2, ncol = 3, dimnames = list(c("clayton", "gumbel"),
                                                               c("0.25", "0.5", "0.75")))
for(i in seq_len(nrow(MultiTable_05))){
  for(j in seq_len(ncol(MultiTable_05))){
    copula = rownames(MultiTable_05)[i]
    tau = as.numeric(colnames(MultiTable_05)[j])
    MultiTable_05[i,j] = OneStudyMulti(copula = copula, tau = tau, seed = 100 + i*j,
                                       alpha = 0.05, beta = 0.05)
  }
}

#Vulnerability CoVaR for alpha = beta = 0.05
VulnerTable_05 = matrix(NA, nrow = 2, ncol = 3, dimnames = list(c("clayton", "gumbel"),
                                                                c("0.25", "0.5", "0.75")))
for(i in seq_len(nrow(VulnerTable_05))){
  for(j in seq_len(ncol(VulnerTable_05))){
    copula = rownames(VulnerTable_05)[i]
    tau = as.numeric(colnames(VulnerTable_05)[j])
    VulnerTable_05[i,j] = OneStudyVulnerability(copula = copula, tau = tau, seed = 10000 + i*j,
                                                alpha = 0.05, beta = 0.05)
  }
}

#Aggregating results
fin05 = rbind(as.vector(BivTable_05), as.vector(MultiTable_05), as.vector(VulnerTable_05))
rownames(fin05) = c("Bivariate", "Multi", "Vulnerability")
colnames(fin05) = rep(c("Clayton", "Gumbel"),3)

#Bivariate CoVaR for alpha = beta = 0.01
BivTable_01 = matrix(NA, nrow = 2, ncol = 3, dimnames = list(c("clayton", "gumbel"),
                                                             c("0.25", "0.5", "0.75")))
for(i in seq_len(nrow(BivTable_01))){
  for(j in seq_len(ncol(BivTable_01))){
    copula = rownames(BivTable_01)[i]
    tau = as.numeric(colnames(BivTable_01)[j])
    BivTable_01[i,j] = OneStudyBiv(copula = copula, tau = tau, seed = 10 + i*j,
                                   alpha = 0.01, beta = 0.01)
  }
}

#Multi CoVaR for alpha = beta = 0.01
MultiTable_01 = matrix(NA, nrow = 2, ncol = 3, dimnames = list(c("clayton", "gumbel"),
                                                               c("0.25", "0.5", "0.75")))
for(i in seq_len(nrow(MultiTable_01))){
  for(j in seq_len(ncol(MultiTable_01))){
    copula = rownames(MultiTable_01)[i]
    tau = as.numeric(colnames(MultiTable_01)[j])
    MultiTable_01[i,j] = OneStudyMulti(copula = copula, tau = tau, seed = 1000 + i*j,
                                       alpha = 0.01, beta = 0.01)
  }
}

#Vulnerability CoVaR for alpha = beta = 0.01
VulnerTable_01 = matrix(NA, nrow = 2, ncol = 3, dimnames = list(c("clayton", "gumbel"),
                                                                c("0.25", "0.5", "0.75")))
for(i in seq_len(nrow(VulnerTable_01))){
  for(j in seq_len(ncol(VulnerTable_01))){
    copula = rownames(VulnerTable_01)[i]
    tau = as.numeric(colnames(VulnerTable_01)[j])
    VulnerTable_01[i,j] = OneStudyVulnerability(copula = copula, tau = tau, seed = 100000 + i*j,
                                                alpha = 0.01, beta = 0.01)
  }
}

#Aggregating results
fin01 = rbind(as.vector(BivTable_01), as.vector(MultiTable_01), as.vector(VulnerTable_01))
rownames(fin01) = c("Bivariate", "Multi", "Vulnerability")
colnames(fin01) = rep(c("Clayton", "Gumbel"),3)

#Final table
round(rbind(fin05, fin01),4)
