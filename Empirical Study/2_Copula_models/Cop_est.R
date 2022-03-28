# setup (loads VaR_violations.R, univariate_model_definition.R and functions.R)
source("C:/Users/MWaltz/Desktop/Forschung/CoVaR/Code/Empirical Study/setup.R")


#------------------------- Copula: Bivariate static est ------------------------

Biv_Cop_n = list()
Biv_Cop_t = list()
Biv_Cop_c = list()
Biv_Cop_g = list()

for(comb in combs){
  for(cop in copula_static){
    
    # which CCs
    cryptos = strsplit(comb, split = "-")[[1]]
    CC1 = cryptos[1]
    CC2 = cryptos[2]
    
    # corresponding fits
    fit1 = fits[[CC1]]
    fit2 = fits[[CC2]]
    
    # estimate copula
    result = get_cop_est_AIC(fit1 = fit1, fit2 = fit2, 
                             copula = cop)
    
    # store results
    switch(cop,
           
           "normal" = {
             Biv_Cop_n[[comb]] = result
           },
           
           "t" = {
             Biv_Cop_t[[comb]] = result
           },
           
           "clayton" = {
             Biv_Cop_c[[comb]] = result
           },
           
           "gumbel" = {
             Biv_Cop_g[[comb]] = result
           },
    )
  }
}

# create table
cop_tbl = rbind(sapply(Biv_Cop_n, function(x) x$param),
                sapply(Biv_Cop_n, function(x) x$AIC),
                sapply(Biv_Cop_t, function(x) x$param[1]),
                sapply(Biv_Cop_t, function(x) x$param[2]),
                sapply(Biv_Cop_t, function(x) x$AIC),
                sapply(Biv_Cop_c, function(x) x$param),
                sapply(Biv_Cop_c, function(x) x$AIC),
                sapply(Biv_Cop_g, function(x) x$param),
                sapply(Biv_Cop_g, function(x) x$AIC))

rownames(cop_tbl) = c("theta_n", "AIC_n", 
                      "theta_t", "df_t", "AIC_t",
                      "theta_c", "AIC_c",
                      "theta_g", "AIC_g")


#----------------------- Copula: Bivariate Patton & DCC ------------------------

Biv_Cop_Pat = list()
Biv_Cop_DCC = list()

# compute time-varying CoVaRs
for(comb in combs){
  for(method in copula_TV){
    
    # which CCs
    cryptos = strsplit(comb, split = "-")[[1]]
    CC1 = cryptos[1]
    CC2 = cryptos[2]
    
    # corresponding fits, specs and ts
    fit1 = fits[[CC1]]
    fit2 = fits[[CC2]]
    
    spec1 = specs[[CC1]]
    spec2 = specs[[CC2]]
    
    ts1 = ts[[CC1]]
    ts2 = ts[[CC2]]
    
    # store result
    switch(method,
           
           "Patton" = {
             Biv_Cop_Pat[[comb]] = get_Patton_est_AIC(fit1 = fit1, fit2 = fit2)
            
           },
           
           "DCC" = {
             Biv_Cop_DCC[[comb]] = get_DCC_biv_est_AIC(uni_spec1 = spec1, 
                                                       uni_spec2 = spec2,
                                                       ts1 = ts1, 
                                                       ts2 = ts2)
           }
    )
  }
}

cop_tbl2 = rbind(sapply(Biv_Cop_Pat, function(x) x$param[1]),
                 sapply(Biv_Cop_Pat, function(x) x$param[2]),
                 sapply(Biv_Cop_Pat, function(x) x$param[3]),
                 sapply(Biv_Cop_Pat, function(x) x$param[4]),
                 sapply(Biv_Cop_Pat, function(x) x$AIC),
                 sapply(Biv_Cop_DCC, function(x) x$param[1]),
                 sapply(Biv_Cop_DCC, function(x) x$param[2]),
                 sapply(Biv_Cop_DCC, function(x) x$param[3]),
                 sapply(Biv_Cop_DCC, function(x) x$AIC))

rownames(cop_tbl2) = c("omega_theta", "beta_theta", "c_theta", "v", "AIC_Pat",
                       "a", "b", "v", "AIC")


#------------------ Copula: Multivariate static estimation ---------------------

h1 = get_cop_est_AIC_higherD(fits_CC, copula = "normal")
h2 = get_cop_est_AIC_higherD(fits_CC, copula = "t")
h3 = get_cop_est_AIC_higherD(fits_CC, copula = "clayton")
h4 = get_cop_est_AIC_higherD(fits_CC, copula = "gumbel")


#--------------------- Copula: Multivariate DCC-t-copula -----------------------

data = do.call(cbind, ts_CC)
result = get_DCC_higherD_est_AIC(list_of_specs = specs_CC, data = data)
