
# execute Rel_exc.R
#source("C:/Users/MWaltz/Desktop/Forschung/CoVaR/Code/Empirical Study/Rel_exc.R")

ylim = list("BTC" = c(-1.25, 0.2),
            "ETH" = c(-1.65, 0.35),
            "LTC" = c(-1.2, 0.6),
            "XMR" = c(-1.3, 0.6),
            "XRP" = c(-1.5, 0.6))

points = list("BTC" = c(-0.7, -0.9, -1.1),
              "ETH" = c(-1.1, -1.3, -1.5),
              "LTC" = c(-0.7, -0.9, -1.1),
              "XMR" = c(-0.8, -1.0, -1.2),
              "XRP" = c(-0.9, -1.1, -1.3))

#-------------------------- Plot: Bivariate CoVaR ------------------------------

plot_BivCoVaR = function(ts1, ts2, name1, name2, VaR1, VaR2, CoVaR, main = F,
                         width = 12, height = 7){
  # setup
  par(cex.axis = 1.25, cex.lab = 1.25)
  
  # init pdf
  pdf(file = paste("BivCoVaR_", name1, "_", name2, ".pdf", sep = ""),
      width = width, height = height)
  
  # prep for axis
  year = 2016:2021
  year_sta = paste(year, "01-01", sep = "-")
  year_idx = which(is.element(date_ts, year_sta))
  
  # evaluate VaR2
  VaR1_eval = eval_VaR(ts1, VaR1)
  VaR2_eval = eval_VaR(ts2, VaR2)
  
  # evaluate
  CoVaR_eval = eval_Biv(ts1 = ts1, ts2 = ts2, VaR2 = VaR2, CoVaR = CoVaR)
  
  #create plot
  if(main){
    main = paste("CoVaR of ", name1, "|", name2,  "- Realized violation rate: ", 
                 CoVaR_eval$rate, sep = "")
  }else{
    main = ""
  }

  plot.ts(cbind(ts1, VaR1, CoVaR), plot.type = "single", col = c("blue3", "red3", "green3"),
          main = main, ylim = ylim[[name1]], ylab = "log-returns", xaxt = "n", lwd = 2)
  
  # add axis
  axis(1, at = year_idx, labels = year)
  
  # cases, when ts1 <= VaR1
  points(VaR1_eval$exc, rep(points[[name1]][1], length(VaR1_eval$exc)), pch = 4)
  
  # cases, when condition is fulfilled (ts2 <= VaR2)
  points(VaR2_eval$exc, rep(points[[name1]][2], length(VaR2_eval$exc)), pch = 4)
  
  # exceedences given the fulfilled condition
  points(CoVaR_eval$CoVaR_exc, rep(points[[name1]][3], length(CoVaR_eval$CoVaR_exc)), pch = 4)
  
  # legend and info
  legend = c(paste("log-return of", name1), paste("5% VaR of", name1), "5% CoVaR")
  legend("bottom", horiz = TRUE, legend = legend,
         fill = c("blue3", "red3", "green3"), bty = "n")
  
  text(200, points[[name1]][1] + 0.05, bquote(.(name1)<=VaR(.(name1))~":"), pos = 4, cex = 1.05)
  text(200, points[[name1]][2] + 0.05, bquote(.(name2)<=VaR(.(name2))~":"), pos = 4, cex = 1.05)
  text(200, points[[name1]][3] + 0.05, bquote(.(name1)<=CoVaR~"|"~.(name2)<=VaR(.(name2))~":"), 
       pos = 4, cex = 1.05)
  
  # close file
  dev.off()
}

# get relevant combs
#combs_rel = combs[startsWith(combs, "BTC")]
combs_rel = combs

for(comb in combs_rel){

  # which CCs
  cryptos = strsplit(comb, split = "-")[[1]]
  CC1 = cryptos[1]
  CC2 = cryptos[2]
  
  # do not consider SCoVaR (is in separate plot below)
  if(startsWith(CC2, "Sys")){
    next
  }
  
  # get ts and VaR
  ts1 = ts[[CC1]]
  ts2 = ts[[CC2]]
  
  VaR1 = VaR[[CC1]]
  VaR2 = VaR[[CC2]]
  
  # create plot
  plot_BivCoVaR(ts1 = ts1, ts2 = ts2, name1 = CC1, name2 = CC2, 
                VaR1 = VaR1, VaR2 = VaR2, CoVaR = Biv_CoVaRs_t[[comb]])
}


#------------------------- Plot: SCoVaR|MCoVaR|VCoVaR --------------------------

ylim_VSM = list("BTC" = c(-3.175, 0.2),
                "ETH" = c(-3.175, 0.3),
                "LTC" = c(-3.175, 0.5),
                "XMR" = c(-3.175, 0.6),
                "XRP" = c(-3.175, 0.7))

plot_VSM_CoVaR = function(name_Y, main = F, width = 12, height = 7){
  
  # get information
  ts_Y   = ts[[name_Y]]
  VaR_Y  = VaR[[name_Y]]
  MCoVaR = MCoVaRs_t[[name_Y]]
  VCoVaR = VCoVaRs_t[[name_Y]]
  SCoVaR = Biv_CoVaRs_t[[which(grepl("Sys", names(Biv_CoVaRs_t)) & 
                                 startsWith(names(Biv_CoVaRs_t), "BTC"))]]
  ts_X  = ts_CC[-which(names(ts_CC) == name_Y)]
  VaR_X = VaR_CC[-which(names(VaR_CC) == name_Y)]

  # setup
  par(cex.axis = 1.25, cex.lab = 1.25)
  
  # init pdf
  pdf(file = paste("VSM_CoVaR_", name_Y, ".pdf", sep = ""),
      width = width, height = height)
  
  # prep for axis
  year = 2016:2021
  year_sta = paste(year, "01-01", sep = "-")
  year_idx = which(is.element(date_ts, year_sta))
  
  # evaluate VaR_Y
  VaR_Y_eval = eval_VaR(ts_Y, VaR_Y)
  
  # evaluate VCoVaR, SCoVaR, MCoVaR
  VCoVaR_eval = eval_MCoVaR_VCoVaR(measure = "VCoVaR", ts_Y = ts_Y, 
                                   MCoVaR_VCoVaR = VCoVaR, ts_X = ts_X, 
                                   VaR_X = VaR_X)
  
  SCoVaR_eval = eval_Biv(ts1   = ts_Y, 
                         ts2   = eval(parse(text = paste("rSys", name_Y, sep = ""))),
                         CoVaR = SCoVaR, 
                         VaR2  = eval(parse(text = paste("VaR_rSys", name_Y, sep = ""))))
  
  MCoVaR_eval = eval_MCoVaR_VCoVaR(measure = "MCoVaR", ts_Y = ts_Y, 
                                   MCoVaR_VCoVaR = MCoVaR, ts_X = ts_X, 
                                   VaR_X = VaR_X)
  
  # plot/text/point setup
  if(main){
    main = paste("VCoVaR, SCoVaR, MCoVaR of", name_Y)
  }else{
    main = ""
  }
  
  col   = c("blue3", "red3", rgb(0, 0.675, 0), "darkgreen", "skyblue4")
  at    = -1.125 - cumsum(c(0, rep(c(0.08, 0.28), 6)))
  coord = 1500
  pos   = 2
  
  # create plot
  plot.ts(cbind(ts_Y, VaR_Y, VCoVaR, SCoVaR, MCoVaR), plot.type = "single", 
          col = col, main = main, ylim = ylim_VSM[[name_Y]], ylab = "log-returns", 
          xaxt = "n", lwd = 1.5)
 
  # add axis
  axis(1, at = year_idx, labels = year)
  
  # VCoVaR condition
  text(coord, at[1], expression(symbol("\044")~i:~ X[i]<="VaR("~X[i]~"):"), 
       pos = pos, cex = 1)
  
  points(VCoVaR_eval$exc_X_agg, rep(at[2], length(VCoVaR_eval$exc_X_agg)), 
         pch = 4)
  
  # VCoVaR exceedences
  text(coord, at[3], bquote(.(name_Y)<=VCoVaR ~"|"~symbol("\044")~i:~ X[i]<="VaR("~X[i]~"):"), 
       pos = pos, cex = 1, col = col[3])

  points(VCoVaR_eval$MCoVaR_VCoVaR_exc, 
         rep(at[4], length(VCoVaR_eval$MCoVaR_VCoVaR_exc)), 
         pch = 4, col = col[3])
  
  # SCoVaR condition
  text(coord, at[5], expression(sum(X[i])<="VaR("~sum(X[i])~"):"), pos = pos, 
       cex = 1)

  points(SCoVaR_eval$exc, rep(at[6], length(SCoVaR_eval$exc)), pch = 4)
  
  # SCoVaR exceedences
  text(coord, at[7], bquote(.(name_Y)<=SCoVaR ~"|"~sum(X[i])<="VaR("~sum(X[i])~"):"), 
       pos = pos, cex = 1, col = col[4])
  
  points(SCoVaR_eval$CoVaR_exc, rep(at[8], length(SCoVaR_eval$CoVaR_exc)), 
         pch = 4, col = col[4])
  
  # MCoVaR condition
  text(coord, at[9], expression(symbol("\042")~i:~ X[i]<="VaR("~X[i]~"):"), 
       pos = pos, cex = 1)
  
  points(MCoVaR_eval$exc_X_agg, rep(at[10], length(MCoVaR_eval$exc_X_agg)), 
         pch = 4)
  
  # MCoVaR exceedences
  text(coord, at[11], bquote(.(name_Y)<=MCoVaR ~"|"~symbol("\042")~i:~ X[i]<="VaR("~X[i]~"):"), 
       pos = pos, cex = 1, col = col[5])

  points(MCoVaR_eval$MCoVaR_VCoVaR_exc, 
         rep(at[12], length(MCoVaR_eval$MCoVaR_VCoVaR_exc)), 
         pch = 4, col = col[5])
  
  # legend
  legend = c(paste("log-return of", name_Y), paste("5% VaR of", name_Y), 
             "5% VCoVaR", "5% SCoVaR", "5% MCoVaR")
  legend("bottom", horiz = TRUE, legend = legend, fill = col, bty = "n")
  
  # close file
  dev.off()
}

plot_VSM_CoVaR(name_Y = "BTC")
plot_VSM_CoVaR(name_Y = "ETH")
plot_VSM_CoVaR(name_Y = "LTC")
plot_VSM_CoVaR(name_Y = "XMR")
plot_VSM_CoVaR(name_Y = "XRP")
