# setup.R and function_oos.R
source("C:/Users/MWaltz/Desktop/Forschung/CoVaR/Code/Empirical Study/setup.R")
source("C:/Users/MWaltz/Desktop/Forschung/CoVaR/Code/Empirical Study/functions_oos.R")

window.length = 500

#--------------------- VaR: OOS forecasting using rolling window ---------------

OOS_VaRs = list()

# compute OOS VaRs
for(name in names(ts)){
    
    # forecast VaR
    OOS_VaRs[[name]] = loop_over_data(data = as.matrix(ts[[name]], ncol = 1), measure = "VaR", 
                                      Y_name = NULL, copula = NULL, window.length = window.length, 
                                      n.ahead = 1, alpha = 0.05, beta = NULL)
}

#---------------------------- VaR: OOS evaluation ------------------------------

OOS_VaRs_eval = list()

for(name in names(ts)){
  OOS_VaRs_eval[[name]] = eval_VaR(ts = ts[[name]][-c(1:window.length)],
                                   VaR = OOS_VaRs[[name]])$rate
}


#-------------------- CoVaR: OOS forecasting using rolling window --------------

Biv_OOS_CoVaRs_n = list()
Biv_OOS_CoVaRs_t = list()
Biv_OOS_CoVaRs_c = list()
Biv_OOS_CoVaRs_g = list()

# compute OOS CoVaRs
for(comb in combs){
  for(cop in copula_static){

    # which CCs
    cryptos = strsplit(comb, split = "-")[[1]]
    CC1 = cryptos[1]
    CC2 = cryptos[2]
    
    # setup data
    data = cbind(ts[[CC1]], ts[[CC2]])
    colnames(data) = cryptos

    # forecast CoVaR
    CoVaR = loop_over_data(data = data, measure = "CoVaR", Y_name = NULL,
                           copula = cop, window.length = window.length, 
                           n.ahead = 1, alpha = 0.05, beta = 0.05)

    # store result
    switch(cop,
           
           "normal" = {
             Biv_OOS_CoVaRs_n[[comb]] = CoVaR
           },
           
           "t" = {
             Biv_OOS_CoVaRs_t[[comb]] = CoVaR
           },
           
           "clayton" = {
             Biv_OOS_CoVaRs_c[[comb]] = CoVaR
           },
           
           "gumbel" = {
             Biv_OOS_CoVaRs_g[[comb]] = CoVaR
           },
    )
  }
}


#------------------------ OOS CoVaR: Bivariate Evaluation ----------------------

# prepare output matrix
out = matrix(nrow = length(combs), ncol = length(copula_static))
rownames(out) = combs
colnames(out) = copula_static

# evaluate
for(i in seq_along(combs)){
  for(j in seq_along(copula_static)) {
    
    # which CCs
    cryptos = strsplit(combs[i], split = "-")[[1]]
    CC1 = cryptos[1]
    CC2 = cryptos[2]
    
    # get ts and OOS-VaR of X
    ts1  = ts[[CC1]][-c(1:window.length)]
    ts2  = ts[[CC2]][-c(1:window.length)]
    VaR2 = OOS_VaRs[[CC2]]
    
    # get CoVaR
    switch(copula_static[j],
           
           "normal" = {
             CoVaR = Biv_OOS_CoVaRs_n[[combs[i]]]
           },
           
           "t" = {
             CoVaR = Biv_OOS_CoVaRs_t[[combs[i]]]
           },

           "clayton" = {
             CoVaR = Biv_OOS_CoVaRs_c[[combs[i]]]
           },
           
           "gumbel" = {
             CoVaR = Biv_OOS_CoVaRs_g[[combs[i]]]
           }
    )
    
    # evaluate
    out[i, j] = eval_Biv(ts1 = ts1, ts2 = ts2, CoVaR = CoVaR, VaR2 = VaR2)$rate
    
  }
}


#----------------- MCoVaR/VCoVaR: OOS forecasting using rolling window ---------

OOS_MCoVaRs_n = list()
OOS_MCoVaRs_t = list()
OOS_MCoVaRs_c = list()
OOS_MCoVaRs_g = list()

OOS_VCoVaRs_n = list()
OOS_VCoVaRs_t = list()
OOS_VCoVaRs_c = list()
OOS_VCoVaRs_g = list()


for(name in names(ts_CC)){
  for(cop in copula_static){
    
    # forecast MCoVaR
    MCoVaR = loop_over_data(data = do.call(cbind, ts_CC), measure = "MCoVaR", 
                            Y_name = name, copula = cop, window.length = window.length, 
                            n.ahead = 1, alpha = 0.05, beta = 0.05)
    
    # forecast VCoVaR
    VCoVaR = loop_over_data(data = do.call(cbind, ts_CC), measure = "VCoVaR", 
                            Y_name = name, copula = cop, window.length = window.length, 
                            n.ahead = 1, alpha = 0.05, beta = 0.05)
    
    # store result
    switch(cop,
           
           "normal" = {
             OOS_MCoVaRs_n[[name]] = MCoVaR
             OOS_VCoVaRs_n[[name]] = VCoVaR
           },
           
           "t" = {
             OOS_MCoVaRs_t[[name]] = MCoVaR
             OOS_VCoVaRs_t[[name]] = VCoVaR
           },
           
           "clayton" = {
             OOS_MCoVaRs_c[[name]] = MCoVaR
             OOS_VCoVaRs_c[[name]] = VCoVaR
           },
           
           "gumbel" = {
             OOS_MCoVaRs_g[[name]] = MCoVaR
             OOS_VCoVaRs_g[[name]] = VCoVaR
           }
    )
  }
}

#------------------------- OOS: MCoVaR/VCoVaR Evaluation -----------------------

# get reduced list of OOS VaRs
OOS_VaRs_CC = OOS_VaRs[names(ts_CC)]

# prepare output matrix
out_MV = matrix(nrow = 2*length(ts_CC), ncol = length(copula_static))
rownames(out_MV) = c(paste(names(ts_CC), c("MCoVaR"), sep = "-"), 
                     paste(names(ts_CC), c("VCoVaR"), sep = "-"))
colnames(out_MV) = copula_static

# evaluate
for(i in seq_len(nrow(out_MV))){
  for(j in seq_along(copula_static)) {
    
    # which crypto and measure
    info = strsplit(rownames(out_MV)[i], split = "-")[[1]]
    CC      = info[1]
    measure = info[2]
    
    # get ts_Y
    ts_Y = ts_CC[[CC]][-c(1:window.length)]
    
    # get list of ts and list of OOS VaRs (X)
    ts_X  = ts_CC[-which(names(ts_CC) == CC)]
    VaR_X = OOS_VaRs_CC[-which(names(OOS_VaRs_CC) == CC)]
    
    # cut first window.length observations of ts_X entries
    for(name in names(ts_X)){
      ts_X[[name]] = ts_X[[name]][-c(1:window.length)]
    }
    
    # get MCoVaR/VCoVaR
    switch(copula_static[j],
           
           "normal" = {
             
             if(measure == "MCoVaR"){
               MCoVaR_VCoVaR = OOS_MCoVaRs_n[[CC]]
             }else{
               MCoVaR_VCoVaR = OOS_VCoVaRs_n[[CC]]
             }
             
           },
           
           "t" = {
             
             if(measure == "MCoVaR"){
               MCoVaR_VCoVaR = OOS_MCoVaRs_t[[CC]]
             }else{
               MCoVaR_VCoVaR = OOS_VCoVaRs_t[[CC]]
             }
             
           },
           
           "clayton" = {
             
             if(measure == "MCoVaR"){
               MCoVaR_VCoVaR = OOS_MCoVaRs_c[[CC]]
             }else{
               MCoVaR_VCoVaR = OOS_VCoVaRs_c[[CC]]
             }
             
           },
           
           "gumbel" = {
             
             if(measure == "MCoVaR"){
               MCoVaR_VCoVaR = OOS_MCoVaRs_g[[CC]]
             }else{
               MCoVaR_VCoVaR = OOS_VCoVaRs_g[[CC]]
             }
             
           }
    )
    
    # evaluate
    out_MV[i, j] = eval_MCoVaR_VCoVaR(measure = measure, ts_Y = ts_Y,
                                      MCoVaR_VCoVaR = MCoVaR_VCoVaR,
                                      ts_X = ts_X, VaR_X = VaR_X)$rate
  }
}
