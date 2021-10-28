Iteratesims_F <- function(df, 
                        Fin = NA,
                        model = 'linear'
                        ){
  
  nruns <- length(Fin)
  
  df.save <- data.frame(years = rep(df$years, nruns),
                        F0 = rep(Fin, each =length(df$years)),
                        SSB = NA,
                        R = NA,
                        Rtot = NA,
                        Catch = NA, 
                        run = rep(1:nruns, each = length(df$years)),
                        model = model)
  
  
  for(i in 1:length(Fin)){

    df$F0 <- df$F0*0+Fin[i] # Just make sure they're the same length
    
    tmprun <- run.agebased.true.catch(df, seed = seeds[i])
    
    
    df.save[df.save$run == i,]$SSB <- tmprun$SSB
    df.save[df.save$run == i,]$R <- tmprun$R.save
    df.save[df.save$run == i,]$Rtot <- tmprun$Rtot.save
    df.save[df.save$run == i,]$Catch <- tmprun$Catch
    
    
    
    
  }
  

  
  
  return(list(ts = df.save))
  
  
  
  
  
}