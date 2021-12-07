runScenarios <- function(models = c('linear','hyper'),
                         recLambda = c('noBOFF','BOFF'),
                         nruns = 100,
                         years = 100,
                         lambda.in = NA,
                         runLambda = TRUE,
                         egg.df,
                         lambda.slope = 0.3,
                         Linf = 150,
                         maxage = 10,
                         tau = 5,
                         K = 0.4,
                         M = 0.4,
                         SDR = .5,
                         F0 = 0,
                         R0 = 1000,
                         h = 0.4,
                         recruitment = 'BH_R',
                         recruitment.type = 'AR',
                         fishing.type = 'constant',
                         mortality = 'constant',
                         seeds = round(runif(years, min = 1, max = 1e6))){
  
  
  
  if(runLambda == TRUE){
    
    
    
    df <- load_data_seasons(nseason = 1,
                            nyear = years,# Set up parameters 
                            Linf = Linf, 
                            maxage = maxage,
                            tau = tau,
                            K = K, 
                            t0 = t0, 
                            M= M,
                            SDR = 0, # Recruitment deviations - set to zero to calculate lambda
                            fishing.type = 'constant',
                            mortality = 'constant',
                            recruitment = recruitment,
                            negg = codest$parameters[['alpha.lin']],
                            eggbeta = codest$parameters[['beta.lin']],
                            F0 = 0, # Set to zero to calc lambda
                            R0 = R0) # Specify parameters
    
    
    tmp <- run.agebased.true.catch(df)
    
    lambda.in <- sum(tmp$N.save.age[df$age >= df$tau,df$tEnd,,]*df$Matsel[df$age >= df$tau]*df$egg.size[df$age >= df$tau])/
      sum(tmp$N.save.age[,df$tEnd,,]*df$Matsel*df$egg.size)
    
    
    print(paste(round(lambda.in,2)*100,'% of spawners is above 50% maturity under no fishing', sep = ''))
    
  }
  
  
  
  
  
  for(k in 1:length(models)){
    
    for(j in 1:length(recLambda)){
      
      df.save <- data.frame(years = rep(1:years, nruns),
                            SSB = NA,
                            R = NA,
                            Rtot = NA,
                            Catch = NA, 
                            run = rep(1:nruns, each = years),
                            model = paste(models[k],recLambda[j], sep = '-'))
      
      
      
      for(i in 1:nruns){
        set.seed(seeds[i])
        
        
        if(models[k] == 'linear'){
          
          negg = egg.df$parameters[['alpha.lin']]
          eggbeta = egg.df$parameters[['beta.lin']]
        }
        
        if(models[k] == 'hyper'){
          
          negg = egg.df$parameters[['alpha.hyper']]
          eggbeta = egg.df$parameters[['beta.hyper']]
        }
        
        if(recLambda[j] == 'noBOFF'){
          lambda = NA
        }
        
        if(recLambda[j] == 'BOFF'){
          lambda = lambda.in
        }
        
        df <- load_data_seasons(nseason = 1,
                                nyear = years,# Set up parameters 
                                Linf = Linf, 
                                maxage = maxage,
                                K = K, 
                                t0 = t0, 
                                M = M,
                                SDR = SDR, # Recruitment deviations 
                                fishing.type = fishing.type,
                                mortality = mortality,
                                alpha = alpha,
                                beta = beta,
                                recruitment = recruitment,
                                recruitment.type = recruitment.type,
                                negg = negg,
                                eggbeta = eggbeta,
                                F0 = F0,
                                R0 = R0,
                                lambda = lambda,
                                lambda.slope = lambda.slope) # Specify parameters
        
        
        
        tmprun <- run.agebased.true.catch(df, seed = seeds[i])
        
        
        df.save[df.save$run == i,]$SSB <- tmprun$SSB
        df.save[df.save$run == i,]$R <- tmprun$R.save
        df.save[df.save$run == i,]$Rtot <- tmprun$Rtot.save
        df.save[df.save$run == i,]$Catch <- tmprun$Catch
        
      }
      
      
      df.sum <- df.save %>% group_by(years, model) %>% 
        summarise(S = median(SSB),
                  Rec = median(R),
                  C = median(Catch),
                  Rtot = median(Rtot),
                  Smin = quantile(SSB, probs = 0.05),
                  Smax = quantile(SSB, probs = 0.95),
                  Recmin = quantile(R, probs = 0.05),
                  Rmax = quantile(R, probs = 0.95),
                  Cmin = quantile(Catch, probs = 0.05),
                  Cmax = quantile(Catch, probs = 0.95)
        )
      
      if(j == 1 & k == 1){
        df.sum.out <- df.sum
        df.save.out <- df.save
      }else{
        
        df.sum.out <- rbind(df.sum.out, df.sum)
        df.save.out <- rbind(df.save.out,df.save)
      }
      
    }
  }
  
  
  
return(list(df.sum = df.sum.out,
           df.save = df.save.out))  
  
  
}