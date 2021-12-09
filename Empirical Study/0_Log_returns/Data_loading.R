#Note: The calculations in this script take only a couple of seconds on an Intel(R) Xeon(R)
#      Gold 6136 CPU with 3.00 GHz on a 64-Bit Windows Server 2016 with 24 processors.

#load CC data
CC_data <- read.csv("C:/Users/MWaltz/Desktop/Forschung/CoVaR/Code/Empirical Study/0_Log_returns/CC_data_v2.csv")

date_ts = CC_data$date

BTC = CC_data$BTC
ETH = CC_data$ETH
LTC = CC_data$LTC
XMR = CC_data$XMR
XRP = CC_data$XRP

#log-returns
rBTC = diff(log(BTC))
rETH = diff(log(ETH))
rLTC = diff(log(LTC))
rXMR = diff(log(XMR))
rXRP = diff(log(XRP))

# system definitions
rSysBTC = rETH + rLTC + rXMR + rXRP
rSysETH = rBTC + rLTC + rXMR + rXRP
rSysLTC = rBTC + rETH + rXMR + rXRP
rSysXMR = rBTC + rETH + rLTC + rXRP
rSysXRP = rBTC + rETH + rLTC + rXMR

#remove CC data
remove(CC_data)
