library("rugarch")
library("WeightedPortTest")

# setup (loads VaR_violations.R, univariate_model_definition.R and functions.R)
source("C:/Users/MWaltz/Desktop/Forschung/CoVaR/Code/Empirical Study/setup.R")


#------------------------------ helper functions -------------------------------

fill_out = function(out, fit, name){
  
  # get names, coef, se of fit
  names_fit = names(coef(fit))
  coef = round(coef(fit), 4)
  se = setNames(round(fit@fit$se.coef, 4), names(coef(fit)))
  
  # get col idx
  c_idx = which(colnames(out) == name)
  
  for(i in seq_along(names_fit)){
    
    # get row idx
    r_idx = which(rownames(out) == names_fit[i])
    
    # fill coef
    out[r_idx, c_idx] = coef[i]
    
    # fill se
    out[r_idx, c_idx + 1] = se[i]
  }
  
  return(out)
}


AIC_test = function(fit, ts, ARMA_df, GARCH_df){
  
  vec = c("AIC" = infocriteria(fit)[1]*length(ts)
          , 
          "LB" = Box.test(as.numeric(residuals(fit, standardize = T)), lag = 8, 
                          fitdf = ARMA_df, type = c("Ljung-Box"))$p.value,
          
          "LBsq" = Box.test(as.numeric(residuals(fit, standardize = T))^2, lag = 8, 
                            fitdf = GARCH_df, type = c("Ljung-Box"))$p.value,
          
          "WLM" = Weighted.LM.test(as.numeric(residuals(fit)), h.t = as.numeric(sigma(fit))^2, 
                                   weighted = T, lag = 8, type = "correlation", fitdf = GARCH_df)$p.value,
          
          "Sign Bias" = signbias(fit)[1,2],
          "Neg. Sign Bias" = signbias(fit)[2,2],
          "Pos. Sign Bias" = signbias(fit)[3,2],
          "Joint effect"= signbias(fit)[4,2])
  
  return(round(vec, 4))
}


#------------------------------- Coef table: CC --------------------------------

# setup names
max_names = c("mu", "ar1", "ar2", "ma1", "ma2", "omega", "alpha1", "alpha2", "alpha3",
              "alpha4", "alpha5", "gamma1", "gamma2", "gamma3", "gamma4", "gamma5", 
              "beta1", "skew", "shape")
cols = c("BTC", "BTC_se", "ETH", "ETH_se", "LTC", "LTC_se", "XMR", "XMR_se",
         "XRP", "XRP_se")

# create out matrix
uni_out_CC = matrix(nrow = length(max_names), ncol = length(cols),
                 dimnames = list(max_names, cols))

for(name in names(fits_CC)){
  uni_out_CC = fill_out(uni_out_CC, fit = fits_CC[[name]], name = name)
}

# make it export ready for LaTeX
out_tex = data.frame(row.names = max_names)

for(c_idx in seq(1, ncol(uni_out_CC), 2)){
  paste(formatC(uni_out_CC[,1], digits = 4, format = "f"))
  out_tex[, ncol(out_tex) + 1] = paste(formatC(uni_out_CC[, c_idx], digits = 4, format = "f"),
                                       " (", 
                                       formatC(uni_out_CC[, c_idx + 1], digits = 4, format = "f"), 
                                       ")", 
                                       sep = "")
}

colnames(out_tex) = names(fits_CC)
out_tex[out_tex == "   NA (   NA)"] = ""

out_tex_CC = out_tex


#----------------------------- Coef table: Systems -----------------------------

# setup names
max_names = c("mu","omega", "alpha1", "beta1", "skew", "shape")
cols = c("SysBTC", "SysBTC_se", "SysETH", "SysETH_se", "SysLTC", "SysLTC_se", 
         "SysXMR", "SysXMR_se", "SysXRP", "SysXRP_se")

# create out matrix
uni_out_sys = matrix(nrow = length(max_names), ncol = length(cols),
                     dimnames = list(max_names, cols))

for(name in c("SysBTC", "SysETH", "SysLTC", "SysXMR", "SysXRP")){
  uni_out_sys = fill_out(uni_out_sys, fit = fits[[name]], name = name)
}

# make it export ready for LaTeX
out_tex = data.frame(row.names = max_names)

for(c_idx in seq(1, ncol(uni_out_sys), 2)){
  paste(formatC(uni_out_sys[,1], digits = 4, format = "f"))
  out_tex[, ncol(out_tex) + 1] = paste(formatC(uni_out_sys[, c_idx], digits = 4, format = "f"),
                                       " (", 
                                       formatC(uni_out_sys[, c_idx + 1], digits = 4, format = "f"), 
                                       ")", 
                                       sep = "")
}

colnames(out_tex) = c("SysBTC", "SysETH", "SysLTC", "SysXMR", "SysXRP")
out_tex[out_tex == "   NA (   NA)"] = ""

out_tex_Sys = out_tex


#---------------------------- AIC and test tables ------------------------------

cbind("BTC" = AIC_test(fit = fitrBTC, ts = rBTC, ARMA_df = 4, GARCH_df = 5),
      "ETH" = AIC_test(fit = fitrETH, ts = rETH, ARMA_df = 0, GARCH_df = 3),
      "LTC" = AIC_test(fit = fitrLTC, ts = rLTC, ARMA_df = 0, GARCH_df = 2),
      "XMR" = AIC_test(fit = fitrXMR, ts = rXMR, ARMA_df = 0, GARCH_df = 2),
      "XRP" = AIC_test(fit = fitrXRP, ts = rXRP, ARMA_df = 0, GARCH_df = 2))

cbind("SysBTC" = AIC_test(fit = fitrSysBTC, ts = rSysBTC, ARMA_df = 0, GARCH_df = 2),
      "SysETH" = AIC_test(fit = fitrSysETH, ts = rSysETH, ARMA_df = 0, GARCH_df = 2),
      "SysLTC" = AIC_test(fit = fitrSysLTC, ts = rSysLTC, ARMA_df = 0, GARCH_df = 2),
      "SysXMR" = AIC_test(fit = fitrSysXMR, ts = rSysXMR, ARMA_df = 0, GARCH_df = 2),
      "SysXRP" = AIC_test(fit = fitrSysXRP, ts = rSysXRP, ARMA_df = 0, GARCH_df = 2))
