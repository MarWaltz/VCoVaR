
#----------------------------- Margin fitting ----------------------------------

fit_margins = function(ts_matrix){
  
  # define general specification
  spec = ugarchspec(variance.model = list(garchOrder = c(1,1), model = "gjrGARCH"),
                    mean.model = list(armaOrder = c(0,0), include.mean = FALSE), 
                    distribution.model = "sstd")
  
  # init output list
  list_of_fits = list()
  
  # estimate each margin
  for(i in seq_len(ncol(ts_matrix))){
    list_of_fits[[i]] = ugarchfit(spec, data = ts_matrix[, i], solver = "hybrid")
  }
  
  return(list_of_fits)
}


#--------------------------- Univariate forecasting ----------------------------

forecast_margin = function(fit, quantile, n.ahead = n.ahead){
  
  # forecast margin
  fore = ugarchforecast(fitORspec = fit, n.ahead = n.ahead)
  
  # calculate risk measure (VaR, CoVaR, etc.)
  quant_OOS = as.vector(quantile(fore, quantile))
  
  return(quant_OOS)
}


#----------------------------- VaR: OOS forecast -------------------------------

get_VaR_OOS = function(ts_matrix, alpha, n.ahead){

  # fit margins
  fit = fit_margins(ts_matrix = ts_matrix)[[1]]
  
  # perform forecast
  VaR_OOS = forecast_margin(fit = fit, quantile = alpha, n.ahead = n.ahead)
  
  return(VaR_OOS)
}


#------------------------ Bivariate CoVaR: OOS forecast ------------------------

get_CoVaR_OOS = function(ts_matrix, copula, n.ahead, alpha, beta){
  # Note: Calculates OOS CoVaR of ts1 given ts2 is in distress.
  
  # fit margins
  list_of_fits = fit_margins(ts_matrix = ts_matrix)
  
  # estimate copula
  param = get_cop_est_AIC(fit1 = list_of_fits[[1]], fit2 = list_of_fits[[2]],
                          copula = copula)$param

  # define objective function
  switch(copula,
         
         "normal" = {
           MinF = function(v){
             return(pCopula(c(v, alpha), copula = normalCopula(param = param, dispstr = "un")) - alpha*beta)
           }
         },
         
         "t" = {
           MinF = function(v){
             return(pCopula(c(v, alpha), copula = tCopula(param = param[1], df = round(param[2]), dispstr = "un")) - alpha*beta)
           }
         },
         
         "clayton" = {
           MinF = function(v){
             return(pCopula(c(v, alpha), copula = claytonCopula(param = param)) - alpha*beta)
           }
         },
         
         "gumbel" = {
           MinF = function(v){
             return(pCopula(c(v, alpha), copula = gumbelCopula(param = param)) - alpha*beta)
           }
         },
         
  )
  
  # find root (= quantile)
  root = uniroot(f = MinF, lower = 0, upper = 1, tol = .Machine$double.eps*0.0000001, 
                 maxiter = 10000)$root
  
  # compute forecast
  CoVaR = forecast_margin(fit = list_of_fits[[1]], quantile = root, n.ahead = n.ahead)
  
  return(CoVaR)
}


#-------------------------- MCoVaR/VCoVaR: OOS forecast ------------------------

get_MCoVaR_VCoVaR_OOS = function(ts_matrix, Y_name, measure, copula, n.ahead, alpha, beta){

  # check
  if(!is.element(measure, c("MCoVaR", "VCoVaR"))){
    stop("Only 'MCoVaR' and 'VCoVaR' are available in this call.")
  }
  
  # get Y_ind
  Y_ind = which(colnames(ts_matrix) == Y_name)
  
  # get dim
  dim = ncol(ts_matrix)
  
  # fit margins
  list_of_fits = fit_margins(ts_matrix = ts_matrix)
  
  # estimate copula
  param = get_cop_est_AIC_higherD(list_of_fits, copula)$param
  
  # define multivariate copula object
  switch(copula,
         
         "normal" = {
           copMulti = normalCopula(param = param, dim = dim, dispstr = "un")
         },
         
         "t" = {
           copMulti = tCopula(param = param[-length(param)], 
                              df = round(param[length(param)]), 
                              dim = dim, dispstr = "un")
         },
         
         "clayton" = {
           copMulti = claytonCopula(param = param, dim = dim)
         },
         
         "gumbel" = {
           copMulti = gumbelCopula(param = param, dim = dim)
         }
  )
  
  set.seed(123)
  
  if(measure == "VCoVaR"){
    
    # rotate copula
    copVul = rotCopula(copMulti)
    
    # define minimization objective
    MinF = function(v){
      vecV = rep(1-alpha, dim)
      vecV[Y_ind] = 1-v
      
      vecA = rep(1-alpha, dim)
      vecA[Y_ind] = 1
      
      return((v - pCopula(vecA, copula = copVul) + pCopula(vecV, copula = copVul)) - 
               (beta * (1 - pCopula(vecA, copula = copVul))))
    }
  }else{
    
    # define minimization objective
    MinF = function(v){
      vecV = rep(alpha, dim)
      vecV[Y_ind] = v
      
      vecA = rep(alpha, dim)
      vecA[Y_ind] = 1
      
      return(pCopula(vecV, copula = copMulti) - pCopula(vecA, copula = copMulti)*beta)
    }
  }
  
  # find root (= quantile)
  root = uniroot(f = MinF, lower = 0, upper = 1, tol = .Machine$double.eps*0.0000001, 
                 maxiter = 10000)$root
  
  # compute forecast
  MCoVaR_VCoVaR = forecast_margin(fit = list_of_fits[[Y_ind]], quantile = root, n.ahead = n.ahead)
  
  return(MCoVaR_VCoVaR)
}


#------------------------------- forecast loop ---------------------------------

loop_over_data = function(data, measure, Y_name = NULL, copula, window.length = 500, 
                          n.ahead = 1, alpha = 0.05, beta = 0.05){
  # checks
  if(!is.element(measure, c("VaR", "CoVaR", "MCoVaR", "VCoVaR"))){
    stop("Only 'VaR', 'CoVaR', 'MCoVaR', and 'VCoVaR' are available in this call.")
  }
  if(!is.matrix(data)){
    stop("'data' should be a matrix.")
  }
  
  # inits
  res = c()
  nobs = nrow(data)
  stop = FALSE; set.seed(123); i = 1
  
  # main loop
  repeat{
    
    # get current estimation indices
    ind = (1 + n.ahead*(i-1)):(window.length + n.ahead*(i-1))
    
    if((max(ind) + n.ahead) > nobs){
      
      if(n.ahead == 1){
        break
      }else{
        n.ahead = nobs - max(ind)
        stop = TRUE
      }
    }
    
    #cat(paste("Start:", min(ind), "End:", max(ind)))
    #cat("\n")
    
    # select current sample
    curr_data = data[ind, ]
    
    # compute forecast
    switch(measure,
           
           "VaR" = {
             value = get_VaR_OOS(ts_matrix = as.matrix(curr_data, ncol = 1), 
                                 alpha = alpha, n.ahead = n.ahead)
           },
           
           "CoVaR" = {
             value = get_CoVaR_OOS(ts_matrix = curr_data, copula = copula,
                                   n.ahead = n.ahead, alpha = alpha, beta = beta)
           },
           
           "MCoVaR" = {
             value = get_MCoVaR_VCoVaR_OOS(ts_matrix = curr_data, Y_name = Y_name,
                                           copula = copula, measure = measure, n.ahead = n.ahead, 
                                           alpha = alpha, beta = beta)
           },
           
           "VCoVaR" = {
             value = get_MCoVaR_VCoVaR_OOS(ts_matrix = curr_data, Y_name = Y_name,
                                           copula = copula, measure = measure, n.ahead = n.ahead, 
                                           alpha = alpha, beta = beta)
           },
    )
    
    # store value
    res = c(res, value)

    # breaking condition if n.ahead != 1
    if(stop){
      break
    }else{
      i = i + 1
    }
  }
  
  return(res)
}
