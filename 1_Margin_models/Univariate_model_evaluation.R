#Note: The calculations in this script take only a couple of seconds/minutes on an Intel(R) 
#      Xeon(R) Gold 6136 CPU with 3.00 GHz on a 64-Bit Windows Server 2016 with 24 processors.

library("rugarch")
library("WeightedPortTest")

# univariate model definition
source("C:/Users/MWaltz/Desktop/Forschung/CoVaR/Code/Empirical Study/1_Margin_models/Univariate_model_definition.R")

#Parameters
round(coef(fitrBTC), 4)
round(coef(fitrETH), 4)
round(coef(fitrLTC), 4)
round(coef(fitrXMR), 4)
round(coef(fitrXRP), 4)

round(coef(fitrSysBTC), 4)
round(coef(fitrSysETH), 4)
round(coef(fitrSysLTC), 4)
round(coef(fitrSysXMR), 4)
round(coef(fitrSysXRP), 4)

#Standard errors
round(fitrBTC@fit$se.coef, 4)
round(fitrETH@fit$se.coef, 4)
round(fitrLTC@fit$se.coef, 4)
round(fitrXMR@fit$se.coef, 4)
round(fitrXRP@fit$se.coef, 4)

round(fitrSysBTC@fit$se.coef, 4)
round(fitrSysETH@fit$se.coef, 4)
round(fitrSysLTC@fit$se.coef, 4)
round(fitrSysXMR@fit$se.coef, 4)
round(fitrSysXRP@fit$se.coef, 4)


#AIC and test tables
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

cbind("BTC" = AIC_test(fit = fitrBTC, ts = rBTC, ARMA_df = 0, GARCH_df = 6),
      "ETH" = AIC_test(fit = fitrETH, ts = rETH, ARMA_df = 0, GARCH_df = 6),
      "LTC" = AIC_test(fit = fitrLTC, ts = rLTC, ARMA_df = 0, GARCH_df = 2),
      "XMR" = AIC_test(fit = fitrXMR, ts = rXMR, ARMA_df = 0, GARCH_df = 2),
      "XRP" = AIC_test(fit = fitrXRP, ts = rXRP, ARMA_df = 0, GARCH_df = 2))

cbind("SysBTC" = AIC_test(fit = fitrSysBTC, ts = rSysBTC, ARMA_df = 3, GARCH_df = 7),
      "SysETH" = AIC_test(fit = fitrSysETH, ts = rSysETH, ARMA_df = 0, GARCH_df = 2),
      "SysLTC" = AIC_test(fit = fitrSysLTC, ts = rSysLTC, ARMA_df = 1, GARCH_df = 5),
      "SysXMR" = AIC_test(fit = fitrSysXMR, ts = rSysXMR, ARMA_df = 2, GARCH_df = 4),
      "SysXRP" = AIC_test(fit = fitrSysXRP, ts = rSysXRP, ARMA_df = 4, GARCH_df = 3))
