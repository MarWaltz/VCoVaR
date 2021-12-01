#Note: The calculations in this script take only a couple of seconds/minutes on an Intel(R) 
#      Xeon(R) Gold 6136 CPU with 3.00 GHz on a 64-Bit Windows Server 2016 with 24 processors.

library("rugarch")

# load data
source("C:/Users/MWaltz/Desktop/Forschung/CoVaR/Code/Empirical Study/0_Log_returns/Data_loading.R")

set.seed(312)

#----------------------- Return Margin Modeling -----------------------
#BTC
specrBTC = ugarchspec(mean.model = list(armaOrder = c(0,0), include.mean = 1),
                      variance.model = list(model = "sGARCH", garchOrder = c(6,0)),
                      distribution.model = "sstd")
fitrBTC = ugarchfit(spec = specrBTC, data = rBTC, solver = "hybrid")

#ETH
specrETH = ugarchspec(mean.model = list(armaOrder = c(0,0), include.mean = 0),
                      variance.model = list(model = "gjrGARCH", garchOrder = c(6,0)),
                      distribution.model = "sstd")
#specrETH = ugarchspec(mean.model = list(armaOrder = c(0,0), include.mean = 0),
#                      variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
#                      distribution.model = "sstd")
fitrETH = ugarchfit(spec = specrETH, data = rETH, solver = "hybrid")

#LTC
specrLTC = ugarchspec(mean.model = list(armaOrder = c(0,0), include.mean = 0),
                      variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
                      distribution.model = "sstd")
fitrLTC = ugarchfit(spec = specrLTC, data = rLTC, solver = "hybrid")

#XMR
specrXMR = ugarchspec(mean.model = list(armaOrder = c(0,0), include.mean = 1),
                      variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
                      distribution.model = "sstd")
fitrXMR = ugarchfit(spec = specrXMR, data = rXMR, solver = "hybrid")

#XRP
specrXRP = ugarchspec(mean.model = list(armaOrder = c(0,0), include.mean = 0),
                      variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
                      distribution.model = "sstd")
fitrXRP = ugarchfit(spec = specrXRP, data = rXRP, solver = "hybrid")

#----------------------- Systems -----------------------

#SysBTC
specrSysBTC = ugarchspec(mean.model = list(armaOrder = c(1,2), include.mean = 1),
                         variance.model = list(model = "gjrGARCH", garchOrder = c(6,1)),
                         distribution.model = "sstd")
fitrSysBTC = ugarchfit(spec = specrSysBTC, data = rSysBTC, solver = "hybrid")

#SysETH
specrSysETH = ugarchspec(mean.model = list(armaOrder = c(0,0), include.mean = 1),
                         variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
                         distribution.model = "sstd")
fitrSysETH = ugarchfit(spec = specrSysETH, data = rSysETH, solver = "hybrid")

#SysLTC
specrSysLTC = ugarchspec(mean.model = list(armaOrder = c(1,0), include.mean = 0),
                         variance.model = list(model = "gjrGARCH", garchOrder = c(5,0)),
                         distribution.model = "sstd")
#specrSysLTC = ugarchspec(mean.model = list(armaOrder = c(0,0), include.mean = 0),
#                         variance.model = list(model = "gjrGARCH", garchOrder = c(1,1)),
#                         distribution.model = "sstd")
fitrSysLTC = ugarchfit(spec = specrSysLTC, data = rSysLTC, solver = "hybrid")

#SysXMR
specrSysXMR = ugarchspec(mean.model = list(armaOrder = c(1,1), include.mean = 1),
                         variance.model = list(model = "sGARCH", garchOrder = c(4,0)),
                         distribution.model = "sstd")
fitrSysXMR = ugarchfit(spec = specrSysXMR, data = rSysXMR, solver = "hybrid")

#SysXRP
specrSysXRP = ugarchspec(mean.model = list(armaOrder = c(2,2), include.mean = 0),
                         variance.model = list(model = "gjrGARCH", garchOrder = c(3,0)),
                         distribution.model = "sstd")
fitrSysXRP = ugarchfit(spec = specrSysXRP, data = rSysXRP, solver = "hybrid")