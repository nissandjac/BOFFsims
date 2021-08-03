## Load the hake data
# year and age input 
load_data_future <- function(df, 
                             bfuture, 
                             catch,
                             sim.data,
                             mortality = 'constant'){
  #'@description Expand the OM data 
  #'@df data frame with old data 
  #'@bfuture bias adjustment 
  #'@catch catch in new year 
  #'@sim.data list holding the previous years OM 
  
  df.new <- df 
  
  df.new$catch.true = TRUE # For future runs calculate fishing mortality rate 
  df.new$Catch <- catch
  df.new$flag_survey <- c(df$flag_survey,1) # Survey every year? 
  df.new$flag_catch <- c(df$flag_catch,1)
  df.new$years <- c(df$years,max(df$years)+1)
  df.new$nyear <- length(df.new$years)
  df.new$tEnd <- length(df.new$years)*df.new$nseason
  
  df.new$wage_catch <- cbind(df$wage_catch,df$wage_catch[,1])
  df.new$wage_survey <- cbind(df$wage_survey,df$wage_survey[,1])
  df.new$wage_mid <- cbind(df$wage_mid,df$wage_mid[,1])
  df.new$wage_ssb <- cbind(df$wage_ssb,df$wage_ssb[,1])
  df.new$b <- c(df$b,bfuture)
  df.new$flag_sel <- c(df.new$flag_sel,0)
  
  Rdevs <- rnorm(n = 1,mean = 0, sd = exp(df$logSDR))
  #Rdevs <- rep(0, yr.future)
  df.new$parms$Rin <- c(df$parms$Rin,Rdevs)
  
  # if(df.new$catch.true == TRUE){ # Remove old fishing mortality, just in case
  #   rm(df.new$F0)
  # }
  
  # Movement is the same as the year before 
  
  ### Add movement to the new years ####
  move.tmp <- array(0, dim = c(df$nspace,df$nage, df$nseason, df.new$nyear))
  move.tmp[,,,1:df.new$nyear-1] <- df$movemat  
  move.tmp[,,,df$nyear] <- df$movemat[,,,df.new$nyear-1]
  df.new$movemat <- move.tmp
  
  
  # Option for natural mortality 
  if(mortality == 'constant'){
    df.new$M0 <- c(df$M0, df$M0[length(df$M0)])
  }
  
  if(mortality == 'rw'){
    df.new$M0 <- c(df$M0, df$M0[length(df$M0)]*exp(rnorm(n = 1, mean = 0, sd = df$SDM)))
    
    df.new$M0[df$M0 > df$M_limit] <- df$M_limit
    
  }
  
  if(mortality == 'increase'){
    df.new$M0 <- c(df$M0, df$M0[df$nyear]+ df$SDM)
    
    df.new$M0[df.new$M0 > df$M_limit] <- df$M_limit
  }
  
  
  if(mortality == 'decrease'){
    
    df.new$M0 <- c(df$M0, df$M0[df$nyear] - df$SDM)
    
    if(any(df.new$M0) < df$M_limit){
      print('natural mortality out of bounds')
    }
    
    df.new$M0[df.new$M0 < df$M_limit] <- df$M_limit
  }
  
  
  if(mortality == 'AR'){ # Natural mortality autocorrelation 
    
  
      df.new$omega <- c(df$omega,
                        df$rho*df$omega[df.new$nyear-1]+sqrt(1-df$rho^2)*rnorm(1,0,sd = sqrt(df$SDM)))
      df.new$M0<-c(df$M0,df$M0[1]*exp(df.new$omega[df.new$nyear]-0.5*(df$SDM^2)))
      
      
    
  }
  
  
  
  

  return(df.new)
  
  }