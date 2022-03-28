library("rugarch")

# univariate model definition
#source("C:/Users/MWaltz/Desktop/Forschung/CoVaR/Code/Empirical Study/1_Margin_models/Univariate_model_definition.R")

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

create_acf_pacf_plot = function(res, name){
  
  # create file
  pdf(file =  paste(name, "_pacf_acf.pdf", sep = ""), width = 10, height = 6)
  
  # fill it
  par(mfrow = c(2,1), mai = c(0.8, 0.82, 0.3, 0.42))
  acf(res, main = "")
  pacf(res, main = "")
  
  # close
  dev.off()
}


#----------------------------- ACF/PACF PLOTS -----------------------------

create_acf_pacf_plot(resBTC, "BTC")
create_acf_pacf_plot(resETH, "ETH")
create_acf_pacf_plot(resLTC, "LTC")
create_acf_pacf_plot(resXMR, "XMR")
create_acf_pacf_plot(resXRP, "XRP")


create_acf_pacf_plot(resSysBTC, "SysBTC")
create_acf_pacf_plot(resSysETH, "SysETH")
create_acf_pacf_plot(resSysLTC, "SysLTC")
create_acf_pacf_plot(resSysXMR, "SysXMR")
create_acf_pacf_plot(resSysXRP, "SysXRP")
