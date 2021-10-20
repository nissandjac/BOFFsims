Iteratesims_F <- function(df, 
                        Fin = NA,
                        recruitment = 'linear',
                        seeds = round(runif(length(dchange),min = 0,1e6))){
  
  nruns <- length(Fin)
  
  df.save <- data.frame(years = rep(df$years, nruns),
                        F0 = rep(Fin, each =length(df$years)),
                        SSB = NA,
                        R = NA,
                        Rtot = NA,
                        Catch = NA, 
                        run = rep(1:nruns, each = length(df$years)),
                        model = recruitment)
  
  
  for(i in 1:length(dchange)){

    df$F0 <- df$F0*0+dchange[i] # Just make sure they're the same length
    
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
              Cmax = quantile(SSB, probs = 0.95)
    )
  
  
  
  return(list(ts = df.save,
              ts.summed = df.sum))
  
  
  
  
  
}