#Note: The calculations in this script takes only a couple of seconds on an Intel(R) Xeon(R)
#      Gold 6136 CPU with 3.00 GHz on a 64-Bit Windows Server 2016 with 24 processors.

library("rugarch")

# univariate model definition
source("C:/Users/MWaltz/Desktop/Forschung/CoVaR/Code/Empirical Study/1_Margin_models/Univariate_model_definition.R")

# get standardized residuals
resBTC = as.vector(residuals(fitrBTC, standardize = T))
resETH = as.vector(residuals(fitrETH, standardize = T))
resLTC = as.vector(residuals(fitrLTC, standardize = T))
resXMR = as.vector(residuals(fitrXMR, standardize = T))
resXRP = as.vector(residuals(fitrXRP, standardize = T))

# same for Systems
resSysBTC = as.vector(residuals(fitrSysBTC, standardize = T))
resSysETH = as.vector(residuals(fitrSysETH, standardize = T))
resSysLTC = as.vector(residuals(fitrSysLTC, standardize = T))
resSysXMR = as.vector(residuals(fitrSysXMR, standardize = T))
resSysXRP = as.vector(residuals(fitrSysXRP, standardize = T))


#----------------------------- ACF/PACF PLOTS -----------------------------
acf(resBTC, main = "ACF - BTC")
pacf(resBTC, main = "PACF - BTC")

acf(resETH, main = "ACF - ETH")
pacf(resETH, main = "PACF - ETH")

acf(resLTC, main = "ACF - LTC")
pacf(resLTC, main = "PACF - LTC")

acf(resXMR, main = "ACF - XMR")
pacf(resXMR, main = "PACF - XMR")

acf(resXRP, main = "ACF - XRP")
pacf(resXRP, main = "PACF - XRP")


#-------------------------------- System ----------------------------------
acf(resSysBTC, main = "ACF - SysBTC")
pacf(resSysBTC, main = "PACF - SysBTC")

acf(resSysETH, main = "ACF - SysETH")
pacf(resSysETH, main = "PACF - SysETH")

acf(resSysLTC, main = "ACF - SysLTC")
pacf(resSysLTC, main = "PACF - SysLTC")

acf(resSysXMR, main = "ACF - SysXMR")
pacf(resSysXMR, main = "PACF - SysXMR")

acf(resSysXRP, main = "ACF - SysXRP")
pacf(resSysXRP, main = "PACF - SysXRP")
