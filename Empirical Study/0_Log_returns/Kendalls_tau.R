#Note: The calculations in this script take only a couple of seconds on an Intel(R) Xeon(R)
#      Gold 6136 CPU with 3.00 GHz on a 64-Bit Windows Server 2016 with 24 processors.

# load data
source("C:/Users/MWaltz/Desktop/Forschung/CoVaR/Code/Empirical Study/0_Log_returns/Data_loading.R")

#Kendall's tau
round(cor(cbind(rBTC, rETH, rLTC, rXMR, rXRP, 
                rSysBTC, rSysETH, rSysLTC, rSysXMR, rSysXRP), 
          method = "kendall"), 4)
