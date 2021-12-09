#Note: The calculations in this script take only a couple of seconds on an Intel(R) Xeon(R)
#      Gold 6136 CPU with 3.00 GHz on a 64-Bit Windows Server 2016 with 24 processors.

########################################## PLOT DATA #############################################

# load data
source("C:/Users/MWaltz/Desktop/Forschung/CoVaR/Code/Empirical Study/0_Log_returns/Data_loading.R")

MaxMinScale = function(ts){
  return((ts - min(ts))/(max(ts)-min(ts)))
}

col = c("blue3", "cornflowerblue", "green3", "red3", "orange")

# setup axis
#par(mfrow = c(1,2), font.axis = 2, font.lab = 2)
par(mfrow = c(1,2), cex.axis = 1.25, cex.lab = 1.25)

year = 2016:2021
year_sta = paste(year, "01-01", sep = "-")
year_idx = which(is.element(date_ts, year_sta))

# plot prices
plot.ts(cbind(MaxMinScale(BTC), MaxMinScale(ETH), MaxMinScale(LTC), MaxMinScale(XMR), MaxMinScale(XRP)), 
        xaxt = "n", plot.type = "single", col = col, ylab = "Standardized Prices", lwd = 2)
axis(1, at = year_idx, labels = year)
legend("topleft", col = col, bty = "n", legend = c("BTC", "ETH", "LTC", "XMR", "XRP"), 
       border = rep(NA, 4), lty = rep(1, 4), density = rep(0, 4), lwd = 2, cex = 1.25)

# plot returns
plot.ts(cbind(rBTC, rETH, rLTC, rXMR, rXRP), plot.type = "single", col = col, 
        ylab = "log-returns", xaxt = "n", lwd = 2)
axis(1, at = year_idx-1, labels = year)
#legend("topright", col = col, bty = "n", legend = c("BTC", "ETH", "LTC", "XMR", "XRP"), 
#       border = rep(NA, 4), lty = rep(1, 4), density = rep(0, 4), lwd = 2, cex = 1.25)
