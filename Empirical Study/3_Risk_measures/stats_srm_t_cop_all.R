
# execute Rel_exc.R
source("C:/Users/MWaltz/Desktop/Forschung/CoVaR/Code/Empirical Study/Rel_exc.R")


#------------------------------- helper fncs -----------------------------------
stats = function(ts){
  
  ts = as.vector(ts)
  
  return(round(c("Min"    = min(ts),
                 "Mean"   = mean(ts),
                 "Median" = median(ts),
                 "Max"    = max(ts),
                 "Sd"     = sd(ts)), 4))
}

fill_stats = function(out_tbl, name){
  
  # VaR
  out_tbl = rbind(out_tbl, stats(VaR_CC[[name]]))
  rownames(out_tbl)[nrow(out_tbl)] = paste(name, "-VaR", sep = "")
  
  # CoVaRs, SCoVaR
  for(comb in combs[startsWith(combs, name)]){
    out_tbl = rbind(out_tbl, stats(Biv_CoVaRs_t[[comb]]))
    rownames(out_tbl)[nrow(out_tbl)] = comb
  }
  
  # MCoVaR
  out_tbl = rbind(out_tbl, stats(MCoVaRs_t[[name]]))
  rownames(out_tbl)[nrow(out_tbl)] = paste(name, "-MCoVaR", sep = "")
  
  # VCoVaR
  out_tbl = rbind(out_tbl, stats(VCoVaRs_t[[name]]))
  rownames(out_tbl)[nrow(out_tbl)] = paste(name, "-VCoVaR", sep = "")
  
  return(out_tbl)
}


#---------------------------- Calculations -------------------------------------

out_tbl = matrix(ncol = 5, nrow = 0,
                 dimnames = list(c(),c("Min", "Mean", "Median", "Max", "Sd")))

for(name in names(ts_CC)){
  out_tbl = fill_stats(out_tbl, name)
}
