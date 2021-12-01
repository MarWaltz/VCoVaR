library("rugarch")
library("rmgarch")
library("copula")
library("nloptr")

#------------------------------- Value at Risk ---------------------------------

get_VaR = function(fit, alpha = 0.05){
  mu = as.vector(fitted(fit))
  sd = as.vector(sigma(fit))
  
  return(mu + sd * qdist("sstd", alpha, 
                         shape = coef(fit)["shape"], 
                         skew = coef(fit)["skew"]))
}


#------------------------------- VaR Evaluation --------------------------------

eval_VaR = function(ts, VaR){
  exc  = which(ts <= VaR)
  rate = length(exc)/length(ts)
  
  return(list("number_exc" = length(exc), "rate" = rate, "exc" = exc))
}

#------------------------- Copula: Bivariate estimation ------------------------

get_cop_est_AIC = function(fit1, fit2, copula){
  
  # get PIT
  u1 = as.vector(pit(fit1))
  u2 = as.vector(pit(fit2))
  data = cbind(u1, u2)
  
  # est cop
  switch (copula,
          
    "normal" = {
      
      est = try(fitCopula(copula = normalCopula(dispstr = "un"), data = data, 
                          estimate.variance = F, method = "ml"))
      
      if (inherits(est, "try-error")) {
        est = try(fitCopula(copula = normalCopula(dispstr = "un"), data = data, 
                            estimate.variance = F, method = "itau"))
      }
    },
    
    "t" = {
      
      est = try(fitCopula(copula = tCopula(dispstr = "un"), data = data, 
                          estimate.variance = F, method = "ml"))
      
      if (inherits(est, "try-error")) {
        est = try(fitCopula(copula = tCopula(dispstr = "un"), data = data, 
                            estimate.variance = F, method = "itau"))
      }
    },
    
    "clayton" = {
      
      est = try(fitCopula(copula = claytonCopula(), data = data, 
                          estimate.variance = F, method = "ml"))
      
      if (inherits(est, "try-error")) {
        est = try(fitCopula(copula = claytonCopula(), data = data, 
                            estimate.variance = F, method = "itau"))
      }
    },
    
    "gumbel" = {
      
      est = try(fitCopula(copula = gumbelCopula(), data = data, 
                          estimate.variance = F, method = "ml"))
      
      if (inherits(est, "try-error")) {
        est = try(fitCopula(copula = gumbelCopula(), data = data, 
                            estimate.variance = F, method = "itau"))
      }
    },
  )
  
  # get param and compute AIC
  param = est@estimate
  
  if(est@method == "maximum likelihood"){
    AIC = as.numeric(2*(length(coef(fit1)) + length(coef(fit2)) + length(coef(est)))
                     -2*(logLik(est) + fit1@fit$LLH + fit2@fit$LLH))
  }else{
    AIC = NULL
  }
  
  return(list("param" = param, "AIC" = AIC))
}


#---------------------- Copula: Bivariate Patton-t-copula ----------------------

Transform = function(x){
  return(tanh(x/2))
}

UpdatePatton = function(thetaInitial, a, b, c, df, data){
  theta = thetaInitial
  
  for(t in 2:nrow(data)){
    if(t <= 10){
      A = mean(qt(data[1:(t-1),1], df = df)*qt(data[1:(t-1),2], df = df))
    }else{
      A = mean(qt(data[(t-10):(t-1),1], df = df)*qt(data[(t-10):(t-1),2], df = df))
    }
    theta[t] = Transform(a + b * theta[t-1] + c * A)
  }
  return(theta)
}

LogLike = function(dataVec, theta, df){
  return(dCopula(u = dataVec, copula = tCopula(param = theta, df = df, dispstr = "un"), log = TRUE))
}

ObjectiveF = function(param, data){
  a = param[1]
  b = param[2]
  c = param[3]
  df = param[4]
  
  # Initialize parameter vector
  startTheta = fitCopula(copula = tCopula(df = df, df.fixed = TRUE, dispstr = "un"), 
                         data = data, estimate.variance = FALSE)@estimate[1]
  theta = rep(startTheta, nrow(data))
  
  # Update parameter vector according to Patton (2006)
  theta = UpdatePatton(thetaInitial = theta, a = a, b = b, c = c, df = df, data = data)
  
  # Calculate Log-Likelihoods
  LL = apply(cbind(data, theta), 1, function(x){LogLike(dataVec = x[1:2], theta = x[3], df = df)})
  
  return(-sum(LL))
}

get_Patton_est_AIC = function(fit1, fit2){
  
  # get PIT
  u1 = as.vector(pit(fit1))
  u2 = as.vector(pit(fit2))
  data = cbind(u1, u2)
  
  # optimization
  out = nloptr(x0 = c(0.5, 0.5, 0.5, 4), eval_f = ObjectiveF, lb = c(-Inf, -Inf, -Inf, 2),
               ub = c(Inf, Inf, Inf, 200), opts = list(algorithm = "NLOPT_LN_COBYLA", xtol_rel = 1e-04),
               data = data)
  
  # resulting estimate
  param = out$solution
  
  # AIC
  AIC_1 = 2*(length(coef(fit1)) + length(coef(fit2)) + length(param))
  AIC_2 = -2*(-out$objective + fit1@fit$LLH + fit2@fit$LLH)
  AIC = AIC_1 + AIC_2
  
  return(list("param" = param, "AIC" = AIC, "data" = data))
}


#---------------------- Copula: Bivariate DCC-t-copula -------------------------

get_DCC_biv_est_AIC = function(uni_spec1, uni_spec2, ts1, ts2){

  # spec
  BivDCC_spec = cgarchspec(uspec = multispec(list(uni_spec1, uni_spec2)), 
                           dccOrder = c(1, 1),
                           distribution.model = list(copula = "mvt", 
                                                     method = "ML", 
                                                     time.varying = TRUE, 
                                                     transformation = "parametric"))

  # fit
  BivDCC_fit = try(cgarchfit(spec = BivDCC_spec, data = cbind(ts1, ts2), 
                             solver = c("hybrid", "solnp")))
  if(inherits(BivDCC_fit, "try-error")){
    return(list("param" = rep(NA, 3), "AIC" = NA))
  }
  
  # copula estimates and AIC
  est = coef(BivDCC_fit, type = "dcc")
  AIC = 2*length(coef(BivDCC_fit)) -2*(likelihood(BivDCC_fit))
  
  return(list("param" = est, "AIC" = AIC, "DCC_fit" = BivDCC_fit))
}


#------------------ Copula: Multivariate estimation ----------------------------

get_cop_est_AIC_higherD = function(list_of_fits, copula){
  
  # get dim
  dim = length(list_of_fits)
  
  # get PIT
  for(i in 1:dim) {
    assign(paste("u", i, sep = ""), as.vector(pit(list_of_fits[[i]])))
  }
  data = eval(parse(text=paste("cbind(", paste0(paste("u", 1:dim, sep = ""), 
                                                collapse = ", "), ")", sep = "")))
  
  # est cop
  switch (copula,
          
          "normal" = {
            
            est = try(fitCopula(copula = normalCopula(dim=dim, dispstr = "un"), 
                                data = data, estimate.variance = F, method = "ml"))
            
            if (inherits(est, "try-error")) {
              est = try(fitCopula(copula = normalCopula(dim=dim, dispstr = "un"), 
                                  data = data, estimate.variance = F, method = "itau"))
            }
          },
          
          "t" = {
            
            est = try(fitCopula(copula = tCopula(dim=dim, dispstr = "un"), 
                                data = data, estimate.variance = F, method = "ml"))
            
            if (inherits(est, "try-error")) {
              est = try(fitCopula(copula = tCopula(dim=dim, dispstr = "un"), 
                                  data = data, estimate.variance = F, method = "itau"))
            }
          },
          
          "clayton" = {
            
            est = try(fitCopula(copula = claytonCopula(dim=dim), data = data, 
                                estimate.variance = F, method = "ml"))
            
            if (inherits(est, "try-error")) {
              est = try(fitCopula(copula = claytonCopula(dim=dim), data = data, 
                                  estimate.variance = F, method = "itau"))
            }
          },
          
          "gumbel" = {
            
            est = try(fitCopula(copula = gumbelCopula(dim=dim), data = data, 
                                estimate.variance = F, method = "ml"))
            
            if (inherits(est, "try-error")) {
              est = try(fitCopula(copula = gumbelCopula(dim=dim), data = data, 
                                  estimate.variance = F, method = "itau"))
            }
          },
  )
  
  # get param and compute AIC
  param = est@estimate
  
  if(est@method == "maximum likelihood"){
    
    AIC_1 = as.numeric(2*(sum(sapply(list_of_fits, function(x)length(coef(x)))) + length(coef(est))))
    AIC_2 = as.numeric(-2*(logLik(est) + sum(sapply(list_of_fits, function(x) x@fit$LLH))))
    AIC = AIC_1 + AIC_2
    
  }else{
    AIC = NULL
  }
  
  return(list("param" = param, "AIC" = AIC))
}


#--------------------- Copula: Multivariate DCC-t-copula -----------------------

get_DCC_higherD_est_AIC = function(list_of_specs, data){
  
  # spec
  DCC_spec = cgarchspec(uspec = multispec(list_of_specs), 
                        dccOrder = c(1, 1),
                        distribution.model = list(copula = "mvt", 
                                                  method = "ML", 
                                                  time.varying = TRUE, 
                                                  transformation = "parametric"))
  
  # fit
  DCC_fit = cgarchfit(spec = DCC_spec, 
                      data = data,
                      solver = c("hybrid", "solnp"))

  # copula estimates and AIC
  est = coef(DCC_fit, type = "dcc")
  AIC = 2*length(coef(DCC_fit)) -2*(likelihood(DCC_fit))
  
  return(list("param" = est, "AIC" = AIC, "DCC_fit" = DCC_fit))
}


#---------------------- CoVaR: Bivariate estimation ----------------------------

get_BivCoVaR = function(fit1, fit2, copula, alpha = 0.05, beta = 0.05){
  # Note: Computes CoVaR of ts1 given ts2 is in distress.
  
  # estimate copula
  param = get_cop_est_AIC(fit1, fit2, copula)$param
  
  # define objective for root finding
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
  
  # find root
  root = uniroot(f = MinF, lower = 0, upper = 1, tol = .Machine$double.eps*0.0000001, 
                 maxiter = 10000)$root
  
  # inverse margin operator
  mu = as.vector(fitted(fit1))
  sd = as.vector(sigma(fit1))
  
  CoVaR = mu + sd * qdist("sstd", root, 
                          shape = coef(fit1)["shape"], 
                          skew = coef(fit1)["skew"])
  
  return(CoVaR)  
}


#-------------------- CoVaR: Bivariate estimation - TV -------------------------

get_RhoPath_Patton = function(param, data){
  
  # extract params
  a = param[1]
  b = param[2]
  c = param[3]
  df = param[4]
  
  # init rho
  start_rho = fitCopula(copula = tCopula(df = df, df.fixed = TRUE, dispstr = "un"), 
                        data = data, estimate.variance = FALSE)@estimate[1]
  rho = rep(start_rho, nrow(data))
  
  # perform Patton update
  rho = UpdatePatton(thetaInitial = rho, a = a, b = b, c = c, df = df, 
                     data = data)
  
  return(rho)
}


get_BivCoVaR_TV = function(method, fit1, fit2, uni_spec1, uni_spec2, ts1, ts2,
                           alpha = 0.05, beta = 0.05){
  # Note: Computes CoVaR of ts1 given ts2 is in distress.
  
  if(!is.element(method, c("Patton", "DCC"))){
    stop("Only 'Patton' and 'DCC' are available for time-varying CoVaR.")
  }
  
  if(method == "Patton"){
    
    # estimate copula
    tmp = get_Patton_est_AIC(fit1, fit2)
    param = tmp$param
    data  = tmp$data
    
    # get all R_t and df
    rho = get_RhoPath_Patton(param, data)  
    df = round(param[4])
    
    # get margin info
    mu = as.vector(fitted(fit1))
    sd = as.vector(sigma(fit1))
    
    shape = coef(fit1)["shape"]
    skew  = coef(fit1)["skew"]
    
  }else{
    
    # estimate copula
    tmp = get_DCC_biv_est_AIC(uni_spec1, uni_spec2, ts1, ts2)
    #param = tmp$param
    DCC_fit = tmp$DCC_fit
    
    # get all R_t and df
    rho = sapply(DCC_fit@mfit$Rt, function(x){x[1,2]})
    df = round(coef(DCC_fit)["[Joint]mshape"])
    
    # get margin info
    mu = as.vector(fitted(DCC_fit)[,1])
    sd = as.vector(sigma(DCC_fit)[,1])
    
    shape = coef(DCC_fit)["[ts1].shape"]
    skew  = coef(DCC_fit)["[ts1].skew"]
  }
  
  # compute CoVaR numerically
  CoVaR = c()
  
  for(t in seq_along(rho)){
    
    # define object for root finding
    MinF = function(v){
      return(pCopula(c(v, alpha), copula = tCopula(param = rho[t], df = df, dispstr = "un")) - alpha*beta)
    }
    
    # find root
    root = uniroot(f = MinF, lower = 0, upper = 1, tol = .Machine$double.eps*0.0000001, 
                   maxiter = 10000)$root
    
    # invert margin
    CoVaR[t] = mu[t] + sd[t] * qdist("sstd", root, shape = shape, skew = skew)
  }
  
  return(CoVaR)
}


#------------------------ CoVaR: MCoVaR & VCoVaR -------------------------------

get_MCoVaR_VCoVaR = function(measure, list_of_fits, Y_name, copula, 
                             alpha = 0.05, beta = 0.05){
  
  # Note: Computes MCoVaR or VCoVaR of ts Y_name as the Y-variable.

  if(!is.element(measure, c("MCoVaR", "VCoVaR"))){
    stop("Only 'MCoVaR' and 'VCoVaR' are available in this call.")
  }
  
  # get Y_ind
  Y_ind = which(names(list_of_fits) == Y_name)
  
  # get dim
  dim = length(list_of_fits)
  
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
  
  # find root
  root = uniroot(f = MinF, lower = 0, upper = 1, tol = .Machine$double.eps*0.0000001, 
                 maxiter = 10000)$root
  
  # inverse margin operator
  fitY = list_of_fits[[Y_name]]
  
  mu = as.vector(fitted(fitY))
  sd = as.vector(sigma(fitY))
  
  MCoVaR_VCoVaR = mu + sd * qdist("sstd", root, 
                                  shape = coef(fitY)["shape"], 
                                  skew = coef(fitY)["skew"])
  
  return(MCoVaR_VCoVaR)
}


#-------------------- CoVaR: MCoVaR / VCoVaR - TV (DCC) ------------------------

get_MCoVaR_VCoVaR_DCC = function(measure, list_of_specs, data, Y_name, 
                                 alpha = 0.05, beta = 0.05){
  
  # Note: Computes time-varying MCoVaR or VCoVaR of ts Y_name as the Y-variable.

  if(!is.element(measure, c("MCoVaR", "VCoVaR"))){
    stop("Only 'MCoVaR' and 'VCoVaR' are available in this call.")
  }
  
  # get Y_ind
  Y_ind = which(names(list_of_specs) == Y_name)
  
  # get dim
  dim = length(list_of_specs)
  
  # set colnames of 'data' to access rmgarch output later
  colnames(data) = paste("ts", 1:dim, sep = "")
  
  # estimate DCC model
  tmp = get_DCC_higherD_est_AIC(list_of_specs = list_of_specs, data = data)
  #param = tmp$param
  DCC_fit = tmp$DCC_fit
  
  # get all R_t and df
  rho = DCC_fit@mfit$Rt
  df = round(coef(DCC_fit)["[Joint]mshape"])
  
  # get margin info
  mu = as.vector(fitted(DCC_fit)[, Y_ind])
  sd = as.vector(sigma(DCC_fit)[, Y_ind])
  
  shape = coef(DCC_fit)[paste("[ts", Y_ind, "].shape", sep = "")]
  skew  = coef(DCC_fit)[paste("[ts", Y_ind, "].skew", sep = "")]
  
  # compute MCoVaR/VCoVaR numerically
  MCoVaR_VCoVaR = c()
  
  for(t in seq_along(rho)){
    
    set.seed(t)
    
    # define copula object
    copMulti = tCopula(param = P2p(rho[[t]]), dim = dim, df = df, 
                       df.fixed = TRUE, dispstr = "un")
    
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

    # find root
    root = uniroot(f = MinF, lower = 0, upper = 1, tol = .Machine$double.eps*0.0000001, 
                   maxiter = 10000)$root
    
    # invert margin
    MCoVaR_VCoVaR[t] = mu[t] + sd[t] * qdist("sstd", root, shape = shape, 
                                             skew = skew)
  }
  return(MCoVaR_VCoVaR)
}


#------------------------- CoVaR: Bivariate Evaluation -------------------------

eval_Biv = function(ts1, ts2, CoVaR, VaR2){

  # check VaR exceedences of ts2
  exc = which(ts2 <= VaR2)
  
  # consider only CoVaR of those observations 
  red_CoVaR = rep(NA, length(CoVaR))
  red_CoVaR[exc] = CoVaR[exc]
  
  # get CoVaR violations
  CoVaR_exc  = which(ts1 <= red_CoVaR)
  CoVaR_rate = length(CoVaR_exc)/length(exc)
  
  return(list("rate"      = round(CoVaR_rate, 4), 
              "abs"       = length(CoVaR_exc), 
              "CoVaR_exc" = CoVaR_exc,
              "exc"       = exc))
}


#---------------------- CoVaR: MCoVaR/VCoVaR Evaluation ------------------------

eval_MCoVaR_VCoVaR = function(measure, ts_Y, MCoVaR_VCoVaR, ts_X, VaR_X){
  
  if(!is.element(measure, c("MCoVaR", "VCoVaR"))){
    stop("Only 'MCoVaR' and 'VCoVaR' are available in this call.")
  }
  
  if(!(is.list(ts_X) & is.list(VaR_X))){
    stop("ts_X and VaR_X should be lists.")
  }
  
  if(!all(names(ts_X) == names(VaR_X))){
    stop("ts and VaR should match.")
  }
  
  # compute individual VaR exceedences of the X variables
  exc_X = list()
  
  for(name in names(ts_X)){
    exc_X[[name]] = which(ts_X[[name]] <= VaR_X[[name]])
  }
  
  # aggregate them
  if(measure == "MCoVaR"){
    exc_X_agg = sort(Reduce(intersect, exc_X))
  }else{
    exc_X_agg = sort(unique(unlist(exc_X)))
  }
  
  # consider only MCoVaR/VCoVaR of those observations 
  red_MCoVaR_VCoVaR = rep(NA, length(MCoVaR_VCoVaR))
  red_MCoVaR_VCoVaR[exc_X_agg] = MCoVaR_VCoVaR[exc_X_agg]
  
  # get violations
  MCoVaR_VCoVaR_exc  = which(ts_Y <= red_MCoVaR_VCoVaR)
  MCoVaR_VCoVaR_rate = length(MCoVaR_VCoVaR_exc)/length(exc_X_agg)
  
  return(list("rate" = round(MCoVaR_VCoVaR_rate, 4), 
              "MCoVaR_VCoVaR_exc" = MCoVaR_VCoVaR_exc,
              "exc_X_agg" = exc_X_agg))
}
