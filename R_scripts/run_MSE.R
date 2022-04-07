run_MSE <- function(df, 
                    sim.data,
                    tsim = 50,
                    hcr = 'constF',
                    M = 'constant',
                    seed = 123,
                    run_EM = TRUE,
                    dirs = file.path('results'), 
                    ...){
  
  #' @description Run an MSE
  #' @param df - parameters for OM 
  #' @param sim.data - OM historical run 
  #' @param tsim - number of future years to simulate 
  #' @hcr harvest control rule to employ
  #' @M natural mortality variation 
  #' @seed set the seed
  #' @run_EM run the TMB assessment? 
  
RUNTIME <- proc.time()

set.seed(seed)

tstep <- 1:(df$nyear+tsim-1)

V.save <- data.frame(biomass = NA,
                     catch = NA,
                     TAC = NA,
                     survey = NA,
                     Mtrue = NA,
                     Fmsy = NA,
                     MSY = NA,
                     S0 = NA,
                     time = rep(tstep,2),
                     SEmin = NA,
                     SEmax = NA,
                     SDR = NA,
                     model = rep(c('EM','OM'), each = length(tstep))
)

parameter.save <- data.frame(name = c('r','K','SDR'), 
                             value = NA, SEmin = NA, 
                             SEmax = NA)

S0 <- sim.data
K <- max(sim.data$Catch)*10
TACsave <- matrix(NA, tsim)
dfsave <- df

for(time in 1:tsim){
  
  
  #print(paste('timestep = ',time))
  # Run the EM
  if(time > 1){
    Cnew <- exp(rnorm(1, mean = 0, sd = 0))
    Cerror <- c(Cerror,exp(rnorm(1, mean = 0, sd = 0.)))
    
    df.new <- load_data_future(df.new, 
                               bfuture = df$bfuture,
                               catch = c(sim.data.new$Catch,TAC*Cnew),sim.data = sim.data.new, 
                               mortality = M)
    sim.data.new <- run.agebased.true.catch(df.new, seed)
  }else{
    
    df.new <- df
    sim.data.new <- sim.data
    Cerror <- exp(rnorm(length(sim.data.new$Catch), mean = 0, sd = 0.))
    
    
    
  }

  MSY.OM <- getOM_MSY_deter(df.new, dirs = dirs) # OM reference points 
  
  
  if(run_EM == TRUE){
    df.tmb <- list(
      Catchobs = as.numeric(sim.data.new$Catch)*Cerror,
      survey = sim.data.new$survey,
      nyear = df.new$nyear,
      logSDsurv = (df.new$parms$logSDsurv+0.001),
      #logq = log(1),
      logSDF = log(0.05)
      
      
    )
    
    Fparms <- as.numeric(sim.data.new$Fout)
    Fparms[Fparms>0.99] <- 0.9
    Fparms[Fparms == 0] <- 0.0001
    
    parms <- list(
      invlogit_Fzero = logit(Fparms),
      logBpred = log(as.numeric(sim.data.new$V.save)),
      logr= log(0.99),
      logK = log(sim.data.new$SSB_0*2),
      logSDB = log(0.1),
      logq = df$logq,
      logSDcatch = log(0.1),
      logSDsurv = log(0.05),
      # 
      logSDF = log(0.1)
    )
    
    # Phase 1 
    obj <-MakeADFun(df.tmb,parms,DLL="runAssessment", random = c('logBpred','invlogit_Fzero'),
                    map = list('logq' = as.factor(NA),
                               'logSDcatch' = as.factor(NA)), silent = TRUE)
    
    
    if(min(obj$report()$B) < 0.001){
      
      parms <- list(
        invlogit_Fzero = log(Fparms),
        logBpred = log(as.numeric(sim.data.new$V.save)),
        logr= log(0.99),
        logK = log(sim.data.new$SSB_0*5),
        logSDB = log(0.1),
        logq = df$logq,
        logSDcatch = log(0.1),
        logSDsurv = log(0.05),
        # 
        logSDF = log(0.1)
      )
      
      # Phase 1 
      obj <-MakeADFun(df.tmb,parms,DLL="runAssessment", random = c('logBpred','invlogit_Fzero'),
                      map = list('logq' = as.factor(NA),
                                 'logSDcatch' = as.factor(NA)), silent = TRUE)
      
      
    }
    
    
    
    lower <- obj$par-Inf
    upper <- obj$par+Inf
    # 
    # upper[names(upper) == 'logr'] <- log(4)
    # upper[names(upper) == 'logSDR'] <- log(2)
    # upper[names(upper) == 'logK'] <- log(max(df.tmb$Catchobs)*5)
    # 
    # 
    lower[names(lower) == 'logr'] <- log(0.1)
    lower[names(lower) == 'logSDB'] <- log(0.01)
    #lower[names(lower) == 'logSDcatch'] <- log(0.01)
    lower[names(lower) == 'logK'] <- log(max(df.tmb$Catchobs)*2.5)
    # 
    upper[names(upper) == 'logr'] <- log(2)
    upper[names(upper) == 'logSDB'] <- log(0.5)
    upper[names(upper) == 'logK'] <- log(max(df.tmb$Catchobs)*15)
    upper[names(upper) == 'logSDcatch'] <- log(0.3)
    upper[names(upper) == 'logSDF'] <- log(1)
    
    #obj$env$tracepar <- TRUE
    
    opt<-nlminb(obj$par,obj$fn,obj$gr,lower=lower,upper=upper, 
                control = list(iter.max = 4000,
                               eval.max = 4000)) #
    reps <- obj$report()
    rep<-TMB::sdreport(obj)
    
    
    if(opt$iterations == 4000){
      print('not converged')
    }
    # 
    if(opt$convergence != 0){
      warning('not converged')
    }

    # Get the TAC  
    if(hcr == 'constF'){
      parameters = data.frame(F0 = Fin)
      
      
      TAC <- HCR_rule(hcr = hcr,biomass = Bhcr,
                      yr = df$tEnd+time-1,parameters = parameters,
                      EM = run_EM)
    }
  
  if(hcr == 'cfp'){
    parameters <- list(Catch = reps$Catch,
                       Bobs = reps$Bpred,
                       K = exp(opt$par)['logK'],
                       r = exp(opt$par)['logr']
                       )
  }
  
  if(hcr == 'cfp' | hcr == 'Fmsy'){
    parameters <- list(Catch = reps$Catch,
                       Bobs = reps$Bpred,
                       K = exp(opt$par)['logK'],
                       r = exp(opt$par)['logr']
    )
    
    Bhcr <- reps$Bpred
    
    
    TAC <- HCR_rule(hcr = hcr,biomass = Bhcr,
                    yr = df$tEnd+time-1,parameters = parameters,
                    EM = run_EM)
  }
  }else{
    if(hcr == 'constF'){
      parameters = data.frame(F0 = Fin)
    }
    
    if(hcr == 'cfp' | hcr == 'Fmsy'){
      parameters <- list(Catch = sim.data.new$Catch,
                         Bobs = sim.data.new$V.save[,1,1],
                         MSY = MSY.OM$MSY,
                         Fmsy = MSY.OM$Fmsy,
                         Bmsy = MSY.OM$Bmsy
      )
    }
    
    Bhcr <- sim.data.new$V.save[,1,1]
    
    TAC <- HCR_rule(hcr = hcr,biomass = Bhcr,
                    yr = df$tEnd+time-1,parameters = parameters,
                    run_EM = run_EM)
    #   plot(sim.data.new$SSB, xlim = c(0,100), ylim = range(sim.data.new$SSB)*c(0.1,2))
    # }else{
    #   points(sim.data.new$SSB)
    #   lines(sim.data.new$SSB)
    # }
  }
  
  

  
  
  if(is.na(TAC)){
    stop('wrong tac')
  }
  
  
  if(is.nan(TAC)){
    stop('wrong tac')
  }
  
  # 
  if(TAC == 0){
    TAC <- sim.data.new$SSB_0*0.01
  }
    
  
    TACsave[time+1] <- TAC
  
    tix <- df$nyear+time-1
    
    if(run_EM == TRUE){
    Fmsy.EM <- as.numeric(exp(opt$par['logr'])*0.5)
    MSY.EM <- as.numeric(exp(opt$par['logr'])*exp(opt$par)['logK']/4)
    

    }
    
    
        
    if(time == 1){
      V.save$Fmsy[V.save$model == 'OM' & V.save$time%in% 1:tix] <- MSY.OM[[2]]
      V.save$MSY[V.save$model == 'OM' & V.save$time%in% 1:tix] <- MSY.OM[[1]]
      V.save$S0[V.save$model == 'OM' & V.save$time%in% 1:tix] <- calcSSB0(df, df$nyear)
      V.save$TAC[V.save$model == 'OM' & V.save$time%in% 1:tix] <- df$Catch      
      
      if(run_EM == TRUE){
      V.save$Fmsy[V.save$model == 'EM' & V.save$time%in% 1:tix] <-  Fmsy.EM
      V.save$MSY[V.save$model == 'EM' & V.save$time%in% 1:tix] <- MSY.EM
      V.save$S0[V.save$model == 'EM' & V.save$time%in% 1:tix] <-  as.numeric(exp(opt$par)['logK'])
      V.save$SDR[V.save$model == 'EM' & V.save$time%in% 1:tix] <- as.numeric(exp(opt$par)['logSDB'])
      
      
      }
      
      }else{
      V.save$Fmsy[V.save$model == 'OM' & V.save$time == tix] <- MSY.OM[[2]]
      V.save$MSY[V.save$model == 'OM' & V.save$time == tix] <- MSY.OM[[1]]
      V.save$S0[V.save$model == 'OM' & V.save$time == tix] <- calcSSB0(df.new, df.new$nyear)
      V.save$TAC[V.save$model == 'OM' & V.save$time == tix+1] <- TAC
      
      
      if(run_EM == TRUE){
      V.save$Fmsy[V.save$model == 'EM' & V.save$time == tix] <-  Fmsy.EM
      V.save$MSY[V.save$model == 'EM' & V.save$time == tix] <- MSY.EM
      V.save$S0[V.save$model == 'EM' & V.save$time == tix] <-  as.numeric(exp(opt$par)['logK'])
      V.save$SDR[V.save$model == 'EM' & V.save$time == tix] <- as.numeric(exp(opt$par)['logSDB'])
      
      }else{
        
      }
    }
    
  
    # }
    
 # Save last years assessment  + uncertainty
  if(time == tsim){
    
    #plot(sim.data.new$R.save, type = 'l')
    
    if(run_EM == TRUE){
    
    rep<-TMB::sdreport(obj)
    
    #print(rep)
    sdrep <- summary(rep)
    rep.values<-rownames(sdrep)
    
    df.Bpred <- getUncertainty('Bpred', df.tmb, sdrep, rep.values)
    Catch <- getUncertainty('Catch', df.tmb, sdrep, rep.values)
    survey <- getUncertainty('surv_est', df.tmb, sdrep, rep.values)
    df.Bpred$model <- 'Estimated'
    
    # Save time varying estimates 
    V.save$biomass = c(df.Bpred$value,as.numeric(sim.data.new$V.save))
    V.save$survey = c(df.Bpred$value,as.numeric(sim.data.new$V.save))
    V.save$SEmin[1:length(tstep)] <- df.Bpred$min
    V.save$SEmax[1:length(tstep)] <- df.Bpred$max
    V.save$catch <- c(Catch$value,sim.data.new$Catch)
    V.save$Mtrue[V.save$model == 'OM'] <- df.new$M0
       # Save parameter estimation 
    nms <- unlist(strsplit(names(rep$par.fixed),  # Remove log from names 
                           split = 'log'))[seq(2,length(names(rep$par.fixed))*2, by = 2)]
    
    parameter.save <- data.frame(name = nms, 
                                 value = NA, SEmin = NA, 
                                 SEmax = NA,
                                 convergence = NA)
    
    for(j in 1:length(rep$par.fixed)){
      tmp <- getUncertainty(names(rep$par.fixed)[j],df.tmb,sdrep,rep.values, convers = 'log')
      
      parameter.save$value[parameter.save$name == nms[j]] <- tmp$value
      parameter.save$SEmin[parameter.save$name == nms[j]] <- tmp$min
      parameter.save$SEmax[parameter.save$name == nms[j]] <- tmp$max
      
    }
    
    }else{
      V.save$biomass[V.save$model == 'OM'] <- as.numeric(sim.data.new$V.save)
      V.save$survey[V.save$model == 'OM'] <- as.numeric(sim.data.new$V.save)
      V.save$catch[V.save$model == 'OM'] <- sim.data.new$Catch
      V.save$Mtrue[V.save$model == 'OM'] <- df.new$M0
      
    }
  }


 
} # End time loop
print(proc.time() - RUNTIME)


# Plot catch
# ggplot(V.save, aes(x = time, y = catch, color = model))+geom_line()+
#   ggplot(V.save, aes(x = time, y = biomass, color = model))+geom_line()

# Recalculate exploitation rate to fishing mortality rate 
#V.save$Fmsy[V.save$model == 'EM'] <- -log(1-V.save$Fmsy[V.save$model == 'EM'])

return(list(timeseries = V.save,
            parameters = parameter.save,
            sim.data = sim.data.new))
}

