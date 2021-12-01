#Note: The calculations in this script take only a couple of seconds on an Intel(R) Xeon(R)
#      Gold 6136 CPU with 3.00 GHz on a 64-Bit Windows Server 2016 with 24 processors.

library("tseries")
library("timeDate")

# load data
source("C:/Users/MWaltz/Desktop/Forschung/CoVaR/Code/Empirical Study/0_Log_returns/Data_loading.R")

TestPipe = function(vec){
  vec = as.vector(vec)
  
  #Jarque-Bera test: H_0: data is normal distributed
  JB = jarque.bera.test(vec)$p.value
  
  #Ljung-Box test: H_0: no serial correlation
  LB = Box.test(vec, lag = 8, type = c("Ljung-Box"))$p.value
  
  #ADF-test: H_0: ts has unit root
  ADF = adf.test(vec, k = 8)$statistic
  
  #PP-test: H_0: ts has unit root
  PP = PP.test(vec)$statistic
  
  #KPSS-test: H_0: ts is (level) stationary
  KPSS = kpss.test(vec)$statistic
  
  return(round(c("Min" = min(vec),
                 "Mean" = mean(vec),
                 "Median" = median(vec),
                 "Max" = max(vec),
                 "Sd" = sd(vec),
                 "Kurtosis" = kurtosis(vec),
                 "Skewness" = skewness(vec),
                 "JB" = JB,
                 "LB" = LB,
                 "ADF" = ADF,
                 "PP" = PP,
                 "KPSS" = KPSS),4))
}

# compute statistics and tests
res = cbind(TestPipe(rBTC), TestPipe(rETH), TestPipe(rLTC), TestPipe(rXMR), TestPipe(rXRP))
colnames(res) = c("BTC", "ETH", "LTC", "XMR", "XRP")
res

# same for systems
res_sys = cbind(TestPipe(rSysBTC), TestPipe(rSysETH), TestPipe(rSysLTC), 
                TestPipe(rSysXMR), TestPipe(rSysXRP))
colnames(res_sys) = c("SysBTC", "SysETH", "SysLTC", "SysXMR", "SysXRP")
res_sys
