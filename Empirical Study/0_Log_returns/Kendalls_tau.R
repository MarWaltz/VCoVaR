# load data
source("C:/Users/MWaltz/Desktop/Forschung/CoVaR/Code/Empirical Study/0_Log_returns/Data_loading.R")

#Kendall's tau
round(cor(cbind(rBTC, rETH, rLTC, rXMR, rXRP, 
                rSysBTC, rSysETH, rSysLTC, rSysXMR, rSysXRP), 
          method = "kendall"), 4)
