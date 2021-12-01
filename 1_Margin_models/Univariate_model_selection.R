#Note: The calculations in this script take approximately 7h on an Intel(R) Xeon(R)
#      Gold 6136 CPU with 3.00 GHz on a 64-Bit Windows Server 2016 with 24 processors.

library("rugarch")
library("WeightedPortTest")

# load data
source("C:/Users/MWaltz/Desktop/Forschung/CoVaR/Code/Empirical Study/0_Log_returns/Data_loading.R")


getAllARMA = function(ret, ar.max = 6, ma.max = 6){
  table = expand.grid(0:ar.max, 0:ma.max, 0:1)
  extraCols = 3
  for(i in seq_len(extraCols)){
    table = cbind(table, rep(NA, nrow(table)))  
  }
  colnames(table) = c("AR", "MA", "Mean", "AIC", "LB", "LBsq")
  
  for(i in seq_len(nrow(table))){
    #spec and fit
    spec = arfimaspec(mean.model = list(armaOrder = c(table$AR[i], table$MA[i]), include.mean = table$Mean[i]),
                      distribution.model = "sstd")
    
    fit = try(arfimafit(spec, data = as.vector(ret), solver = "hybrid"), silent = T)
    if(class(fit) == "try-error"){next}
    
    #AIC
    table$AIC[i] = infocriteria(fit)[1]
    
    #tests
    res = as.vector(residuals(fit))
    df = sum(c(table$AR[i], table$MA[i]))
    
    tmp = try(Box.test(res, lag = max(df+1, 8), fitdf = df, type = c("Ljung-Box"))$p.value, silent = T)
    if(class(tmp)!= "try-error"){
      table$LB[i] = tmp  
    }
    
    tmp = try(Box.test(res^2, lag = max(df+1, 8), fitdf = df, type = c("Ljung-Box"))$p.value, silent = T)
    if(class(tmp)!= "try-error"){
      table$LBsq[i] = tmp  
    }
    
    res = df = NULL
    cat(paste(i, "/", nrow(table), "\n"))
  }
  return(table)
}


getAllARMAGARCH = function(ret, ar.max = 6, ma.max = 6, arch.max = 6, 
                           garch.max = 1, model = "sGARCH"){

  table = expand.grid(0:ar.max, 0:ma.max, 0:1, 0:arch.max, 0:garch.max)
  
  if(model == "sGARCH"){
    extraCols = 3
  }else{
    extraCols = 7
  }
  for(i in seq_len(extraCols)){
    table = cbind(table, rep(NA, nrow(table)))  
  }
  
  if(model == "sGARCH"){
    colnames(table) = c("AR", "MA", "Mean", "ARCH", "GARCH", "AIC", "LB", "WLM")
  }else{
    colnames(table) = c("AR", "MA", "Mean", "ARCH", "GARCH", "AIC", "LB", "WLM",
                        "SB", "NSB", "PSB", "JE")
  }
  
  for(i in seq_len(nrow(table))){
    
    #spec and fit
    if(model == "sGARCH"){
      spec = ugarchspec(mean.model = list(armaOrder = c(table$AR[i], table$MA[i]),include.mean = table$Mean[i]),
                        variance.model = list(model = "sGARCH", garchOrder = c(table$ARCH[i], table$GARCH[i])),
                        distribution.model = "sstd")
    }
    if(model == "gjrGARCH"){
      spec = ugarchspec(mean.model = list(armaOrder = c(table$AR[i], table$MA[i]),include.mean = table$Mean[i]),
                        variance.model = list(model = "gjrGARCH", garchOrder = c(table$ARCH[i], table$GARCH[i])),
                        distribution.model = "sstd")     
    }
    
    fit = try(ugarchfit(spec, data = as.vector(ret), solver = "hybrid"), silent = T)
    if(class(fit) == "try-error"){next}
    
    #AIC
    table$AIC[i] = infocriteria(fit)[1]
    
    #tests
    res = as.numeric(residuals(fit))
    stdres = as.numeric(residuals(fit, standardize = TRUE))
    h.t = as.numeric(sigma(fit))^2
    df = sum(c(table$AR[i], table$MA[i]))
    gdf = sum(c(table$ARCH[i], table$GARCH[i]))
    
    tmp = try(Box.test(stdres, lag = max(df+1, 8), fitdf = df, type = c("Ljung-Box"))$p.value, silent = T)
    if(class(tmp)!= "try-error"){
      table$LB[i] = tmp  
    }
    
    tmp = try(Weighted.LM.test(res, h.t = h.t, weighted = T, lag = max(gdf+1, 8), type = "correlation", fitdf = gdf)$p.value, silent = T)
    if(class(tmp)!= "try-error"){
      table$WLM[i] = tmp  
    }
    
    if(model != "sGARCH"){
      SB = try(signbias(fit)[,2], silent = T)
      if(class(SB)!= "try-error"){
        table$SB[i]  = SB[1]
        table$NSB[i] = SB[2]
        table$PSB[i] = SB[3]
        table$JE[i]  = SB[4]
      }
    }
    
    res = stdres = h.t = df = gdf = NULL
    cat(paste(i, "/", nrow(table), "\n"))
  }
  return(table)
}


select_model = function(ts){
  
  # Step 1: Calculate all ARMA(p,q) models and consider those fulfilling the LB
  ARMA = getAllARMA(ts, ar.max = 6, ma.max = 6)
  ARMA_red = ARMA[which(ARMA$LB > 0.05),]
  
  
  # Step 2: Select best model according to AIC and test for ARCH effects
  if(nrow(ARMA_red) != 0){
    
    best_idx = order(ARMA_red$AIC)[1]
    ARMA_OPT = ARMA_red[best_idx,]
    
    if(ARMA_OPT$LBsq > 0.05){
      cat("Final model is ARMA.\n")
      return(ARMA_OPT)
    }
  }
  
  cat("ARMA spec not sufficient.\n")
  
  
  # Step 3: Calculate all ARMA(p,q)-GARCH(P,Q) models and consider those fulfilling LB, WLM
  ARMA_GARCH = getAllARMAGARCH(ts, ar.max = 6, ma.max = 6, arch.max = 6, garch.max = 1, 
                               model = "sGARCH")
  ARMA_GARCH_red = ARMA_GARCH[which((ARMA_GARCH$LB > 0.05) & (ARMA_GARCH$WLM > 0.05)),]
  
  
  # Step 4: Select best model according to AIC and test for leverage effects
  if(nrow(ARMA_GARCH_red) != 0){
    
    best_idx = order(ARMA_GARCH_red$AIC)[1]
    ARMA_GARCH_OPT = ARMA_GARCH_red[best_idx,]
    
    spec = ugarchspec(mean.model = list(armaOrder = as.numeric(ARMA_GARCH_OPT[1:2]),
                                        include.mean = as.numeric(ARMA_GARCH_OPT[3])),
                      variance.model = list(model = "sGARCH", 
                                            garchOrder = as.numeric(ARMA_GARCH_OPT[4:5])),
                      distribution.model = "sstd")
    fit = ugarchfit(spec = spec, data = ts, solver = "hybrid")
    
    if(all(signbias(fit)[,2] > 0.05)){
      cat("Final model is ARMA-GARCH.\n")
      return(ARMA_GARCH_OPT)
    }
  }
  
  cat("ARMA-GARCH spec not sufficient.\n")
  
  # Steps 5: Calculate all possible ARMA(p,q)-GJR-GARCH(P,Q) models and
  #          consider those fulfilling LB, WLM, SB
  ARMA_GJR_GARCH = getAllARMAGARCH(ts, ar.max = 6, ma.max = 6, arch.max = 6, 
                                   garch.max = 1, model = "gjrGARCH")
  idx = apply(ARMA_GJR_GARCH, 1, function(x){all(x[7:12] > 0.05)})
  ARMA_GJR_GARCH_red = ARMA_GJR_GARCH[idx,]
  
  
  # Step 6: Select best model according to AIC
  if(nrow(ARMA_GJR_GARCH_red) != 0){
    
    best_idx = order(ARMA_GJR_GARCH_red$AIC)[1]
    ARMA_GJR_GARCH_OPT = ARMA_GJR_GARCH_red[best_idx,]
    
    cat("Final model is ARMA-GJR-GARCH.")
    return(ARMA_GJR_GARCH_OPT)
  }
  
  cat("Model selection not successful.")
}

# new data
modBTC = select_model(rBTC)
modETH = select_model(rETH)
modLTC = select_model(rLTC)
modXMR = select_model(rXMR)
modXRP = select_model(rXRP)

# Optimal models:
# BTC: ARCH(6) with mu
# ETH: GJR-ARCH(6) without mu
# LTC: GARCH(1,1) without mu
# XMR: GARCH(1,1) with mu
# XRP: GARCH(1,1) without mu

modSysBTC = select_model(rSysBTC)
modSysETH = select_model(rSysETH)
modSysLTC = select_model(rSysLTC)
modSysXMR = select_model(rSysXMR)
modSysXRP = select_model(rSysXRP)

# Optimal models:
# SysBTC: ARMA(1,1)-GJR-GARCH(6,1) with mu
# SysETH: GARCH(1,1) with mu
# SysLTC: AR(1)-GJR-GARCH(5,0) without mu
# SysXMR: ARMA(1,1)-ARCH(4,0) with mu
# SysXRP: ARMA(2,2)-GJR-ARCH(3,0) without mu
