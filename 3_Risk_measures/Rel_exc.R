
# setup (loads VaR_violations.R, univariate_model_definition.R and functions.R)
source("C:/Users/MWaltz/Desktop/Forschung/CoVaR/Code/Empirical Study/setup.R")

#------------------------------ Bivariate CoVaR --------------------------------

Biv_CoVaRs_n = list()
Biv_CoVaRs_t = list()
Biv_CoVaRs_c = list()
Biv_CoVaRs_g = list()

# compute CoVaRs
for(comb in combs){
  for(cop in copula_static){

    # which CCs
    cryptos = strsplit(comb, split = "-")[[1]]
    CC1 = cryptos[1]
    CC2 = cryptos[2]
    
    # corresponding fits
    fit1 = fits[[CC1]]
    fit2 = fits[[CC2]]
    
    # estimate CoVaR
    CoVaR = get_BivCoVaR(fit1, fit2, copula = cop)
    
    # store result
    switch(cop,
           
           "normal" = {
             Biv_CoVaRs_n[[comb]] = CoVaR
           },
           
           "t" = {
             Biv_CoVaRs_t[[comb]] = CoVaR
           },
           
           "clayton" = {
             Biv_CoVaRs_c[[comb]] = CoVaR
           },
           
           "gumbel" = {
             Biv_CoVaRs_g[[comb]] = CoVaR
           },
    )
  }
}


#----------------------- CoVaR: Bivariate Patton & DCC -------------------------

Biv_CoVaRs_Pat = list()
Biv_CoVaRs_DCC = list()

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
    
    # estimate CoVaR
    CoVaR = get_BivCoVaR_TV(method = method, fit1 = fit1, fit2 = fit2, 
                            uni_spec1 = spec1, uni_spec2 = spec2, 
                            ts1 = ts1, ts2 = ts2)
    
    # store result
    switch(method,
           
           "Patton" = {
             Biv_CoVaRs_Pat[[comb]] = CoVaR
           },
           
           "DCC" = {
             Biv_CoVaRs_DCC[[comb]] = CoVaR
           }
    )
  }
}

#------------------------- CoVaR: Bivariate Evaluation -------------------------

# prepare output matrix
out = matrix(nrow = length(combs), ncol = length(copula_full))
rownames(out) = combs
colnames(out) = copula_full

# evaluate
for(i in seq_along(combs)){
  for(j in seq_along(copula_full)) {
    
    # which CCs
    cryptos = strsplit(combs[i], split = "-")[[1]]
    CC1 = cryptos[1]
    CC2 = cryptos[2]
    
    # get ts and VaR of X
    ts1  = ts[[CC1]]
    ts2  = ts[[CC2]]
    VaR2 = VaR[[CC2]]
    
    # get CoVaR
    switch(copula_full[j],
           
           "normal" = {
             CoVaR = Biv_CoVaRs_n[[combs[i]]]
           },
           
           "t" = {
             CoVaR = Biv_CoVaRs_t[[combs[i]]]
           },
           
           "Patton" = {
             CoVaR = Biv_CoVaRs_Pat[[combs[i]]]
           },
           
           "DCC" = {
             CoVaR = Biv_CoVaRs_DCC[[combs[i]]]
           },
           
           "clayton" = {
             CoVaR = Biv_CoVaRs_c[[combs[i]]]
           },

           "gumbel" = {
             CoVaR = Biv_CoVaRs_g[[combs[i]]]
           }
    )
    
    # evaluate
    out[i, j] = eval_Biv(ts1 = ts1, ts2 = ts2, CoVaR = CoVaR, VaR2 = VaR2)$rate
    
  }
}


#--------------------------- CoVaR: MCoVaR/VCoVaR ------------------------------

MCoVaRs_n = list()
MCoVaRs_t = list()
MCoVaRs_c = list()
MCoVaRs_g = list()

VCoVaRs_n = list()
VCoVaRs_t = list()
VCoVaRs_c = list()
VCoVaRs_g = list()

for(name in names(fits_CC)){
  for(cop in copula_static){
    
    #if(cop != "t"){next}

    # estimate MCoVaR
    MCoVaR = get_MCoVaR_VCoVaR(measure = "MCoVaR", list_of_fits = fits_CC, 
                               Y_name = name, copula = cop)

    # estimate VCoVaR
    VCoVaR = get_MCoVaR_VCoVaR(measure = "VCoVaR", list_of_fits = fits_CC,
                               Y_name = name, copula = cop)
    
    # store result
    switch(cop,
           
           "normal" = {
             MCoVaRs_n[[name]] = MCoVaR
             VCoVaRs_n[[name]] = VCoVaR
           },
           
           "t" = {
             MCoVaRs_t[[name]] = MCoVaR
             VCoVaRs_t[[name]] = VCoVaR
           },
           
           "clayton" = {
             MCoVaRs_c[[name]] = MCoVaR
             VCoVaRs_c[[name]] = VCoVaR
           },
           
           "gumbel" = {
             MCoVaRs_g[[name]] = MCoVaR
             VCoVaRs_g[[name]] = VCoVaR
           }
    )
  }
}


#--------------------------- CoVaR: DCC MCoVaR/VCoVaR --------------------------

data = do.call(cbind, ts_CC)

MCoVaRs_DCC = list()
VCoVaRs_DCC = list()

for(name in names(specs_CC)){
  
  # estimate MCoVaR
  MCoVaR = get_MCoVaR_VCoVaR_DCC(measure = "MCoVaR", list_of_specs = specs_CC,
                                 data = data, Y_name = name)
    
  # estimate VCoVaR
  VCoVaR = get_MCoVaR_VCoVaR_DCC(measure = "VCoVaR", list_of_specs = specs_CC,
                                 data = data, Y_name = name)
    
  # store result
  MCoVaRs_DCC[[name]] = MCoVaR
  VCoVaRs_DCC[[name]] = VCoVaR
}


#------------------------ CoVaR: MCoVaR/VCoVaR Evaluation ----------------------

# prepare output matrix
out_MV = matrix(nrow = 2*length(ts_CC), ncol = length(copula_full))
rownames(out_MV) = c(paste(names(ts_CC), c("MCoVaR"), sep = "-"), 
                     paste(names(ts_CC), c("VCoVaR"), sep = "-"))
colnames(out_MV) = copula_full

# evaluate
for(i in seq_len(nrow(out_MV))){
  for(j in seq_along(copula_full)) {
    
    # which crypto and measure
    info = strsplit(rownames(out_MV)[i], split = "-")[[1]]
    CC      = info[1]
    measure = info[2]
    
    # get ts_Y
    ts_Y = ts_CC[[CC]]
    
    # get list of ts and list of VaR (X)
    ts_X  = ts_CC[-which(names(ts_CC) == CC)]
    VaR_X = VaR_CC[-which(names(VaR_CC) == CC)]
    
    # get MCoVaR/VCoVaR
    switch(copula_full[j],
           
           "normal" = {
             
             if(measure == "MCoVaR"){
               MCoVaR_VCoVaR = MCoVaRs_n[[CC]]
             }else{
               MCoVaR_VCoVaR = VCoVaRs_n[[CC]]
             }
             
           },
           
           "t" = {
             
             if(measure == "MCoVaR"){
               MCoVaR_VCoVaR = MCoVaRs_t[[CC]]
             }else{
               MCoVaR_VCoVaR = VCoVaRs_t[[CC]]
             }
             
           },
           
           "clayton" = {
             
             if(measure == "MCoVaR"){
               MCoVaR_VCoVaR = MCoVaRs_c[[CC]]
             }else{
               MCoVaR_VCoVaR = VCoVaRs_c[[CC]]
             }
             
           },
           
           "gumbel" = {
             
             if(measure == "MCoVaR"){
               MCoVaR_VCoVaR = MCoVaRs_g[[CC]]
             }else{
               MCoVaR_VCoVaR = VCoVaRs_g[[CC]]
             }
             
           },
           
           "Patton" = {
             
             next
             
           },
           
           "DCC" = {
             
             if(measure == "MCoVaR"){
               MCoVaR_VCoVaR = MCoVaRs_DCC[[CC]]
             }else{
               MCoVaR_VCoVaR = VCoVaRs_DCC[[CC]]
             }
             
           }
    )
    
    # evaluate
    out_MV[i, j] = eval_MCoVaR_VCoVaR(measure = measure, ts_Y = ts_Y,
                                      MCoVaR_VCoVaR = MCoVaR_VCoVaR,
                                      ts_X = ts_X, VaR_X = VaR_X)$rate
  }
}
