library("rugarch")

# univariate model definition
source("C:/Users/MWaltz/Desktop/Forschung/CoVaR/Code/Empirical Study/1_Margin_models/Univariate_model_definition.R")

# functions
source("C:/Users/MWaltz/Desktop/Forschung/CoVaR/Code/Empirical Study/functions.R")


#------------------------------- Value at Risk ---------------------------------
VaR_rBTC = get_VaR(fitrBTC)
VaR_rETH = get_VaR(fitrETH)
VaR_rLTC = get_VaR(fitrLTC)
VaR_rXMR = get_VaR(fitrXMR)
VaR_rXRP = get_VaR(fitrXRP)

VaR_rSysBTC = get_VaR(fitrSysBTC)
VaR_rSysETH = get_VaR(fitrSysETH)
VaR_rSysLTC = get_VaR(fitrSysLTC)
VaR_rSysXMR = get_VaR(fitrSysXMR)
VaR_rSysXRP = get_VaR(fitrSysXRP)

#test via internal VaR calculation of 'rugarch': 
# all(VaR_rBTC == quantile(fitrBTC, 0.05)) #TRUE


#------------------------------- VaR Evaluation --------------------------------

eval_BTC_VaR = eval_VaR(rBTC, VaR_rBTC)
eval_ETH_VaR = eval_VaR(rETH, VaR_rETH)
eval_LTC_VaR = eval_VaR(rLTC, VaR_rLTC)
eval_XMR_VaR = eval_VaR(rXMR, VaR_rXMR)
eval_XRP_VaR = eval_VaR(rXRP, VaR_rXRP)

eval_SysBTC_VaR = eval_VaR(rSysBTC, VaR_rSysBTC)
eval_SysETH_VaR = eval_VaR(rSysETH, VaR_rSysETH)
eval_SysLTC_VaR = eval_VaR(rSysLTC, VaR_rSysLTC)
eval_SysXMR_VaR = eval_VaR(rSysXMR, VaR_rSysXMR)
eval_SysXRP_VaR = eval_VaR(rSysXRP, VaR_rSysXRP)

#---------------------------- aggregate results --------------------------------

res = cbind("alpha-level" = rep(0.05, 10),
            "Measured Rate" = round(c(eval_BTC_VaR[["rate"]],
                                      eval_ETH_VaR[["rate"]],
                                      eval_LTC_VaR[["rate"]],
                                      eval_XMR_VaR[["rate"]],
                                      eval_XRP_VaR[["rate"]],
                                      eval_SysBTC_VaR[["rate"]],
                                      eval_SysETH_VaR[["rate"]],
                                      eval_SysLTC_VaR[["rate"]],
                                      eval_SysXMR_VaR[["rate"]],
                                      eval_SysXRP_VaR[["rate"]]),4),
            "Theoretical Violations" = rep(length(rBTC), 10) * 0.05,
            "Measured Violations" = round(c(eval_BTC_VaR[["number_exc"]],
                                            eval_ETH_VaR[["number_exc"]],
                                            eval_LTC_VaR[["number_exc"]],
                                            eval_XMR_VaR[["number_exc"]],
                                            eval_XRP_VaR[["number_exc"]],
                                            eval_SysBTC_VaR[["number_exc"]],
                                            eval_SysETH_VaR[["number_exc"]],
                                            eval_SysLTC_VaR[["number_exc"]],
                                            eval_SysXMR_VaR[["number_exc"]],
                                            eval_SysXRP_VaR[["number_exc"]]),4))
rownames(res) = c("BTC", "ETH", "LTC", "XMR", "XRP",
                  "SysBTC", "SysETH", "SysLTC", "SysXMR", "SysXRP")
