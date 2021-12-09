
# execute Rel_exc.R
source("C:/Users/MWaltz/Desktop/Forschung/CoVaR/Code/Empirical Study/Rel_exc.R")

# helper fnc
stats = function(ts){
  
  ts = as.vector(ts)
  
  return(round(c("Min"    = min(ts),
                 "Mean"   = mean(ts),
                 "Median" = median(ts),
                 "Max"    = max(ts),
                 "Sd"     = sd(ts)), 4))
}


# get relevant combs
combs_rel = combs[startsWith(combs, "BTC")]

# create/fill out object
out_tbl = stats(VaR_rBTC)
for(comb in combs_rel){
  out_tbl = cbind(out_tbl, stats(Biv_CoVaRs_t[[comb]]))
}

# add MCoVaR/VCoVaR
out_tbl = cbind(out_tbl, 
                stats(MCoVaRs_t[["BTC"]]), 
                stats(VCoVaRs_t[["BTC"]]))

colnames(out_tbl) = c("VaR(BTC)", 
                      paste("CoVaR|", 
                            unlist(strsplit(combs_rel, "-"))[seq(2, length(combs_rel) * 2, 2)],
                            sep = ""),
                      "MCoVaR",
                      "VCoVaR")
