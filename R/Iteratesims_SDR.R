Iteratesims_SDR <- function(df, 
                            SDRin = NA,
                            model = 'linear'
){
  
  nruns <- length(SDRin)
  
  df.save <- data.frame(years = rep(df$years, nruns),
                        SDR = rep(SDRin, each =length(df$years)),
                        SSB = NA,
                        R = NA,
                        Rtot = NA,
                        Catch = NA, 
                        run = rep(1:nruns, each = length(df$years)),
                        model = model)
  
  
  for(i in 1:length(SDRin)){
    
    df$logSDR<- log(SDRin[i])
    
    tmprun <- run.agebased.true.catch(df, seed = seeds[i])
    
    
    df.save[df.save$run == i,]$SSB <- tmprun$SSB
    df.save[df.save$run == i,]$R <- tmprun$R.save
    df.save[df.save$run == i,]$Rtot <- tmprun$Rtot.save
    df.save[df.save$run == i,]$Catch <- tmprun$Catch
    
    
    
    
  }
  
  return(list(ts = df.save))
  
}