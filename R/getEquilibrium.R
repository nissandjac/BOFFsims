getEquilibrium <- function(df){ # 
  
  Msel <- df$Msel # no difference between males and females
  M0 <- df$M0
  
  M <- t(replicate(Msel, n = df$nyear))
  M <- M[1:df$nyear,]*M0
  
  
  # Age 
  nage <- df$nage
  age <- df$age
  
  R0 <- exp(df$parms$logRinit)
  
  Mage <- cumsum(M[1,])#c(0,cumsum(M[1,1:(nage-1)]))
  
  # Calculate N0 based on R0
  mage <- nage # Max age
  agetmp <- 0:(mage)
  nagetmp <- mage*1# SS multiplies this by 3 
  Mtmp <- rep(NA, nagetmp)
  
  N0 <- rep(NA,nagetmp)
  
  N0[1] <- R0
  N0[2:(nagetmp-1)] = R0*exp(-Mage[1:(nagetmp-2)])
  N0[nagetmp] =  R0*exp(-Mage[nagetmp-1])/(1-exp(-M[1,nage]))
  
  # SSB_0 <- sum(N0*df$wage_ssb[,1]*df$Matsel)
  
  # Calculate the R0 required to match the egg production # 
  
  R0new <- sum(N0*df$Matsel*df$egg.size)
  
  return(R0new)
}