# Run egg trial thingie # 
fn_sims <- function(models = c('noBOFF','BOFF'),
                    recLambda = c('hyper','linear'),
                    years = 100,
                    nruns = 100,
                    tau = 5,
                    Linf = 150,
                    maxage = 10,
                    K = 0.6,
                    t0 = 0,
                    SDR = 0.5,
                    SDF = 0,
                    F0 = 0.2,
                    M = 0.4,
                    R0 = 1000,
                    rho = 0.001,
                    egg.scale = 1,
                    tau_sel = 2,
                    Fpast = 0,
                    recruitment = 'BH_steep',
                    lambda.cut = .7,
                    lambda.in = 0.4,
                    mortality = 'constant',
                    fishing.type = 'constant',
                    recruitment.type = 'AR'
                    ){

  
# Just use cod for fun
  eggs <- read.csv('data/fecundityEggSizeFemaleSize.csv')
  eggs$weight <- 0.01*(eggs$FemaleSize_mm/10)^3 # Fix this later
  # All eggs as a function of size 
  
  eggs <- eggs %>% group_by(Species) %>% 
    mutate(relweight = weight/max(weight),                                     
           releggs = Fecundity_nOfEggs_per_female/max(Fecundity_nOfEggs_per_female)
  )
  
  
  x <- eggs[is.na(eggs$relweight) == 0,]$relweight
  y <- eggs[is.na(eggs$relweight) == 0,]$releggs
  
  
  # all relative eggs 
  #parms <- est_eggs(x,y)
  
  
  
  # Just take the cod eggs most of the other fish are not really big fisheries 
  cod <- eggs[eggs$Species == 'Gadus morhua',]
  
  cod$weight <- 0.01*(cod$FemaleSize_mm/10)^3 # Fix the parameters for weight lenght to whatever here 
  
  codest <- est_eggs(x = cod$weight,
                     y = cod$Fecundity_nOfEggs_per_female)
  
  ### Calculate Fmsy 
  # 
  # Fmsy <- seq(0, 3, length.out = 20)
  # 
  # 
  # 
  # for(i in 1:length(Fmsy)){
  # # This df is just used for later calcs 
  # df <- load_data_seasons(nseason = 1,
  #                         nyear = 100,# Set up parameters 
  #                         Linf = Linf, 
  #                         maxage = maxage,
  #                         tau = tau,
  #                         K = K, 
  #                         t0 = t0, 
  #                         M= M,
  #                         tau_sel = tau_sel,
  #                         SDR = 0, # Recruitment deviations - set to zero to calculate lambda
  #                         fishing.type = 'constant',
  #                         mortality = 'constant',
  #                         recruitment = recruitment,
  #                         negg = codest$parameters[['alpha.lin']]/egg.scale,
  #                         eggbeta = codest$parameters[['beta.lin']],
  #                         F0 = Fmsy[i], # Set to zero to calc lambda
  #                         R0 = R0) # Specify parameters
  # 
  # 
  # xx <- run.agebased.true.catch(df)
  # 
  # }
  # 
  
  Fin <- F0
  
  
  ls.plot <- runScenarios(models = models,
                          recLambda = recLambda,
                          nruns = nruns, 
                          years = years,
                          Fpast = Fpast,
                          runLambda = FALSE,
                          lambda.in = lambda.in,
                          rho = rho,
                          egg.df = codest, 
                          egg.scale = egg.scale,
                          lambda.cut = lambda.cut,
                          SDR = SDR,
                          SDF = SDF,
                          F0 = Fin,
                          maxage = maxage,
                          K = K, 
                          Linf = Linf,
                          t0 = t0,
                          tau = tau,
                          tau_sel = tau_sel,
                          M = M,
                          mortality = mortality,
                          fishing.type = fishing.type,
                          recruitment.type = recruitment.type
                          )
  
  # Run all the models 
  

  
  
  # Make Mikaels figures # 
  
  df.N <- ls.plot[[3]]
  df.N <- df.N[df.N$age > 0,]
  
  df.N$SSB <- df.N$N*df.N$weight*df.N$mat
  
  
  ## Summarise the data frame to  plot it 
  df.N$old <- NA
  pold <- tau+2
  
  df.N$old[df.N$age < pold] <- 'young'
  df.N$old[df.N$age >= pold] <- 'old'
  
  # Remove the zero age from the calculations 
  
  
  df.Nsum <- df.N[df.N$age > 0,] %>% 
    group_by(years, model, old, run) %>% 
    dplyr::summarise(N = sum(N),
              Catch = sum(Catch),
              SSB=  sum(SSB),
              F0 = mean(F0),
              Rdev = mean(Rdev)) %>% dplyr::arrange(run, old,model) # The Rdev mean is not a mean (it's the same for all ) 
  
  # Get median weight weighted by numbers
  
  df.wSum <- df.N[df.N$age > 0,] %>% 
    group_by(years, model ,run) %>%
    dplyr::summarise(mWeight = weighted.mean(weight, N),
              mAge = weighted.mean(age, N)) %>% arrange(run, model)
  
  
  
  
  # Rearrange ls.plot[22]

  
  R.df <- ls.plot[[2]] %>% arrange(run, model,years) 
  
  ### DO it the old slow way, but change later 
  modelruns <- unique(R.df$model)
  nruns <- max(R.df$run)
  R.df$residuals <- NA
  
  scammodels <- list()
  
  for(i in 1:nruns){
    for(j in 1:length(modelruns)){
      
      
     dftmp <- R.df[R.df$run == i & R.df$model == modelruns[j],]  
     Ftmp <- scam(R+.001 ~ s(SSB, k = 20, bs = 'mpd', m = 2) +  
                                      offset(log(SSB+.001)), 
                  family=gaussian(link="identity"), data = dftmp,optimizer="nlm",sp=0.01)
     
      
     # Calculate the residuals (anomalies) 
     dftmp$Rpred <- predict(Ftmp, newdata = data.frame(SSB = dftmp$SSB))
     dftmp$residuals <- log(dftmp$Rpred)-log(dftmp$R)
    
     
     # Add the residuals to the R.df data frame 
     R.df[R.df$run == i & R.df$model == modelruns[j],]$residuals <- dftmp$residuals
     
     
     }
  }
  
  
  df.propOld <- data.frame(propOld = df.Nsum$N[df.Nsum$old == 'old']/(df.Nsum$N[df.Nsum$old == 'young']+df.Nsum$N[df.Nsum$old == 'old']),
                           SSBprop = df.Nsum$SSB[df.Nsum$old == 'old']/(df.Nsum$SSB[df.Nsum$old == 'young']+df.Nsum$SSB[df.Nsum$old == 'old']),
                           years = df.Nsum$years[df.Nsum$old == 'old'],
                           model = df.Nsum$model[df.Nsum$old == 'old'],
                           run = df.Nsum$run[df.Nsum$old == 'old'],
                           Rdev = df.Nsum$Rdev[df.Nsum$old == 'old'],
                           Rresid = R.df$residuals,
                           mweight = df.wSum$mWeight,
                           mage = df.wSum$mAge,
                           rec = R.df$R,
                           rho = R.df$rho,
                           lambda = R.df$lambda,
                           M = R.df$M,
                           F0 = R.df$F0,
                           recR0 = R.df$R/exp(df$parms$logRinit),
                           rtot = R.df$Rtot
  )
  
  print(median(df.propOld$SSBprop))
  #ggplot(df.propOld, aes(x = Rresid, y = Rdev))+geom_point()+geom_smooth(method = 'lm')
  
  # 
  # 
  # df.propSum <- df.propOld %>% 
  #   group_by( years, model) %>% 
  #   summarise(propN = median(propOld),
  #             propSSB = median(SSBprop),
  #             mWeight = median(mweight),
  #             mAge = median(mage)) %>% 
  #   pivot_longer(cols = 3:6)
  # 
  # 
  # 
  # ggplot(df.propSum[df.propSum$years > 50,], aes(x = years, y= value, color = model))+geom_line()+facet_wrap(~name, scales = 'free_y')
  
  # Try the correlation coefficients 
  
  prop.plot <- df.propOld %>% group_by(model, run, rho, lambda) %>% dplyr::summarise(SSBcorRR0  = cor(SSBprop, recR0),
                                                                weightcorRR0  = cor(mweight, recR0),
                                                                agecorRR0  = cor(mage, recR0),
                                                                SSBcorRtot  = cor(SSBprop, rtot),
                                                                weightcorRtot  = cor(mweight, rtot),
                                                                agecorRtot  = cor(mage, rtot),
                                                                SSBcorRdev  = cor(SSBprop, Rresid),
                                                                weightcorRdev  = cor(mweight, Rresid),
                                                                agecorRdev  = cor(mage, Rresid)) %>%
    pivot_longer(4:12)
  
  
  lims <- max(abs(prop.plot$value))
 
  
  print(
    ggplot(prop.plot, aes(x = model, y = value, fill = model))+geom_violin()+geom_boxplot(width = 0.1)+
          facet_wrap(~name)+scale_x_discrete()+theme(axis.text.x = element_blank())+
      theme_bw()+coord_cartesian(ylim = c(-lims,lims))+geom_hline(aes(yintercept = 0), linetype = 2)
                                                                        )
  # Create a data frame to save 
  
  df.propOld <- df.propOld %>% mutate(Linf = Linf, SDR = SDR, K = K, tau = tau)
  
  df.export <- df.propOld
  
  return(list(calcs = df.export,
              df.save = ls.plot[[2]],
              N = ls.plot[[3]]))
}



