
# VaR violations (loads univariate_model_definition.R and functions.R)
source("C:/Users/MWaltz/Desktop/Forschung/CoVaR/Code/Empirical Study/1_Margin_models/VaR_violations.R")
rm(res)

#----------------------------------- Setup -------------------------------------

copula_static = c("normal", "t", "clayton", "gumbel")
copula_TV     = c("Patton", "DCC")
copula_full   = c(copula_static, copula_TV)

combs = c("BTC-ETH", "BTC-LTC", "BTC-XMR", "BTC-XRP", "BTC-SysBTC",
          "ETH-BTC", "ETH-LTC", "ETH-XMR", "ETH-XRP", "ETH-SysETH",
          "LTC-BTC", "LTC-ETH", "LTC-XMR", "LTC-XRP", "LTC-SysLTC",
          "XMR-BTC", "XMR-ETH", "XMR-LTC", "XMR-XRP", "XMR-SysXMR",
          "XRP-BTC", "XRP-ETH", "XRP-LTC", "XRP-XMR", "XRP-SysXRP")

fits = list("BTC" = fitrBTC, 
            "ETH" = fitrETH, 
            "LTC" = fitrLTC, 
            "XMR" = fitrXMR, 
            "XRP" = fitrXRP, 
            "SysBTC" = fitrSysBTC, 
            "SysETH" = fitrSysETH,
            "SysLTC" = fitrSysLTC, 
            "SysXMR" = fitrSysXMR, 
            "SysXRP" = fitrSysXRP)

fits_CC = list("BTC" = fitrBTC, 
               "ETH" = fitrETH, 
               "LTC" = fitrLTC, 
               "XMR" = fitrXMR, 
               "XRP" = fitrXRP)

ts = list("BTC" = rBTC, 
          "ETH" = rETH, 
          "LTC" = rLTC, 
          "XMR" = rXMR, 
          "XRP" = rXRP, 
          "SysBTC" = rSysBTC, 
          "SysETH" = rSysETH,
          "SysLTC" = rSysLTC, 
          "SysXMR" = rSysXMR, 
          "SysXRP" = rSysXRP)

ts_CC = list("BTC" = rBTC, 
             "ETH" = rETH, 
             "LTC" = rLTC, 
             "XMR" = rXMR, 
             "XRP" = rXRP)

specs = list("BTC" = specrBTC, 
             "ETH" = specrETH, 
             "LTC" = specrLTC, 
             "XMR" = specrXMR, 
             "XRP" = specrXRP, 
             "SysBTC" = specrSysBTC, 
             "SysETH" = specrSysETH,
             "SysLTC" = specrSysLTC, 
             "SysXMR" = specrSysXMR, 
             "SysXRP" = specrSysXRP)

specs_CC = list("BTC" = specrBTC, 
                "ETH" = specrETH, 
                "LTC" = specrLTC, 
                "XMR" = specrXMR, 
                "XRP" = specrXRP)

VaR = list("BTC" = VaR_rBTC, 
           "ETH" = VaR_rETH, 
           "LTC" = VaR_rLTC, 
           "XMR" = VaR_rXMR, 
           "XRP" = VaR_rXRP, 
           "SysBTC" = VaR_rSysBTC, 
           "SysETH" = VaR_rSysETH,
           "SysLTC" = VaR_rSysLTC, 
           "SysXMR" = VaR_rSysXMR, 
           "SysXRP" = VaR_rSysXRP)

VaR_CC = list("BTC" = VaR_rBTC, 
              "ETH" = VaR_rETH, 
              "LTC" = VaR_rLTC, 
              "XMR" = VaR_rXMR, 
              "XRP" = VaR_rXRP)
