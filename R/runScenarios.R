runScenarios <- function(models = c('linear','hyper'),
                         recLambda = c('noBOFF','BOFF'),
                         nruns = 100,
                         years = 100,
                         lambda.in = NA,
                         runLambda = TRUE,
                         egg.df,
                         lambda.cut = 0.9,
                         Linf = 150,
                         t0 = 0,
                         maxage = 10,
                         tau = 5,
                         K = 0.4,
                         M = 0.4,
                         Fpast = 0,
                         SDR = .5,
                         SDF = 0,
                         F0 = 0,
                         R0 = 1000,
                         h = 0.4,
                         rho = 0.001,
                         egg.scale = 1, 
                         tau_sel= 2,
                         recruitment = 'BH_steep',
                         recruitment.type = 'AR',
                         fishing.type = 'AR',
                         mortality = 'constant',
                         seeds = round(runif(nruns, min = 1, max = 1e6))){
  
  
  
  if(runLambda == TRUE){
    
    
    
    df <- load_data_seasons(nseason = 1,
                            nyear = years,# Set up parameters 
                            Linf = Linf, 
                            maxage = maxage,
                            tau = tau,
                            tau_sel = tau_sel,
                            K = K, 
                            t0 = t0, 
                            M= M,
                            Fpast = Fpast,
                            SDR = 0, # Recruitment deviations - set to zero to calculate lambda
                            fishing.type = fishing.type,
                            mortality = mortality,
                            recruitment = recruitment,
                            recruitment.type = recruitment.type,
                            negg = codest$parameters[['alpha.lin']]/egg.scale,
                            eggbeta = codest$parameters[['beta.lin']],
                            F0 = 0, # Set to zero to calc lambda
                            R0 = R0) # Specify parameters
    
    
    tmp <- run.agebased.true.catch(df)
    
    lambda.in <- sum(tmp$N.save.age[df$age >= df$tau,df$tEnd,,]*df$Matsel[df$age >= df$tau]*df$egg.size[df$age >= df$tau])/
      sum(tmp$N.save.age[,df$tEnd,,]*df$Matsel*df$egg.size)
    
    
    print(paste(round(lambda.in,2)*100,'% of spawners is above 50% maturity under no fishing', sep = ''))
    
  }
  
  
    
 for(k in 1:length(models)){
  
     for(s in 1:length(rho)){
        
       for(l in 1:length(lambda.in)){ 
       
       for(j in 1:length(recLambda)){
      
         
     #    for(p in 1:length(F0)){
           
      
      
         df.save <- data.frame(years = rep(1:years, nruns),
                            F0 = F0,
                            rho = rho[s],
                            SSB = NA,
                            R = NA,
                            Rscam = NA,
                            Rscam_dev = NA,
                            xfrac = NA,
                            R0_boff = NA,
                            SR = NA,
                            Rtot = NA,
                            Catch = NA, 
                            M = NA,
                            lambda = lambda.in[l],
                            Rdev = NA,
                            run = rep(1:nruns, each = years),
                            model = paste(models[k],recLambda[j], sep = '-'))
      
      
        df.N <- data.frame( years = rep(rep(1:years, nruns), each = maxage+1),
                            F0 = NA,
                            rho = NA,
                            lambda = NA,
                            age = rep(0:maxage, nruns*length(1:years)),
                            N = NA,
                            Rdev = NA,
                            R = NA,
                            weight = NA,
                            mat = NA,
                            Catch = NA, 
                            run = rep(rep(1:nruns, each = years), each = maxage+1),
                            model = rep(paste(models[k],recLambda[j], sep = '-'), each = maxage +1)
        )
      
      
      
      
      for(i in 1:nruns){
        
        set.seed(seeds[i])
        
        
        if(models[k] == 'linear'){
          
          negg = egg.df$parameters[['alpha.lin']]/egg.scale
          eggbeta = egg.df$parameters[['beta.lin']]
        }
        
        if(models[k] == 'hyper'){
          
          negg = egg.df$parameters[['alpha.hyper']]/egg.scale
          eggbeta = egg.df$parameters[['beta.hyper']]
        }
        
        if(recLambda[j] == 'noBOFF'){
          lambda = NA
        }
        
        if(recLambda[j] == 'BOFF'){
          lambda = lambda.in[l]
        }
        
        
        if(length(F0) == 1){
          F0 <- rep(F0, years)
        }
        
        
        if(length(F0) != years){
          stop('wrong length of fishing mortality')
        }
        
        df <- load_data_seasons(nseason = 1,
                                nyear = years,# Set up parameters 
                                Linf = Linf, 
                                maxage = maxage,
                                K = K, 
                                t0 = t0, 
                                M = M,
                                tau_sel = tau_sel,
                                tau = tau,
                                SDR = SDR, # Recruitment deviations 
                                fishing.type = fishing.type,
                                mortality = mortality,
                                alpha = alpha,
                                beta = beta,
                                recruitment = recruitment,
                                recruitment.type = recruitment.type,
                                negg = negg,
                                SDF = SDF,
                                eggbeta = eggbeta,
                                F0 = F0,
                                rhoR = rho[s],
                                R0 = R0,
                                lambda = lambda,
                                lambda.cut = lambda.cut) # Specify parameters

        tmprun <- run.agebased.true.catch(df, seed = seeds[i])
        
        
        # Do the scam in here 
        df.tmp <- data.frame(SSB = as.numeric(tmprun$SSB),
                             R = as.numeric(tmprun$R.save))
        
        
        
        Rtmp <- scam(log(R) ~ s(SSB, k = 20, bs = 'mpd', m = 2) +  
                       offset(log(SSB)), 
                     family=gaussian(link="identity"),data = df.tmp, optimizer="nlm",sp=0.01)
        
        
        # Calculate the residuals (anomalies) 
        Rpred <- exp(predict(Rtmp, newdata = data.frame(SSB = df.tmp$SSB)))
        residuals <- log(Rpred)-log(df.tmp$R)
        
        # Calculate fraction of Old and Young
        #young <- tmprun$N.save[*
        
        
        df.save[df.save$run == i,]$SSB <- tmprun$SSB
        df.save[df.save$run == i,]$R <- tmprun$R.save
        df.save[df.save$run == i,]$R0_boff <- tmprun$R0.boff
        df.save[df.save$run == i,]$xfrac <- tmprun$x.save
        df.save[df.save$run == i,]$Rtot <- tmprun$Rtot.save
        df.save[df.save$run == i,]$Catch <- tmprun$Catch
        df.save[df.save$run == i,]$F0 <- df$F0
        df.save[df.save$run == i,]$M <- df$M0
        df.save[df.save$run == i,]$SR <- tmprun$SR
        df.save[df.save$run == i,]$rho <- df$rho
        df.save[df.save$run == i,]$lambda <- df$lambda
        df.save[df.save$run == i,]$Rdev <- df$parms$Rin
        df.save[df.save$run == i,]$Rscam <- Rpred
        df.save[df.save$run == i,]$Rscam_dev <- residuals
        
        df.N[df.N$run == i,]$N <- as.numeric(tmprun$N.save.age[,1:years,,])
        df.N[df.N$run == i,]$weight <- as.numeric(df$wage_ssb)
        df.N[df.N$run == i,]$mat <- rep(df$Matsel, years)
        df.N[df.N$run == i,]$Catch <- as.numeric(tmprun$Catch.age)
        df.N[df.N$run == i,]$rho <- rep(df$rho, each = maxage + 1)
        df.N[df.N$run == i,]$F0 <- as.numeric(tmprun$Fseason)
        df.N[df.N$run == i,]$Rdev <- rep(df$parms$Rin, each = maxage+1)
        df.N[df.N$run == i,]$R <- rep(tmprun$R.save, each = maxage+1)
        df.N[df.N$run == i,]$lambda <- rep(lambda, each = maxage+1)
        
      }
      
      
      df.sum <- df.save %>% group_by(years, model, rho, lambda) %>% 
        summarise(S = median(SSB),
                  Rec = median(R),
                  C = median(Catch),
                  Rtot = median(Rtot),
                  Smin = quantile(SSB, probs = 0.05),
                  Smax = quantile(SSB, probs = 0.95),
                  Recmin = quantile(R, probs = 0.05),
                  Rmax = quantile(R, probs = 0.95),
                  Cmin = quantile(Catch, probs = 0.05),
                  Cmax = quantile(Catch, probs = 0.95),
                  F0 = max(F0)
        )
      
      
      
      if(j == 1 & k == 1 & s == 1 & l == 1){
        df.sum.out <- df.sum
        df.save.out <- df.save
        df.N.out <- df.N
      }else{
        
        df.sum.out <- rbind(df.sum.out, df.sum)
        df.save.out <- rbind(df.save.out,df.save)
        df.N.out <- rbind(df.N.out, df.N)
      }
      
  #  }
       }
       }
  }
  
  }
  
return(list(df.sum = df.sum.out,
           df.save = df.save.out,
           df.N  = df.N.out))  
  
  
}