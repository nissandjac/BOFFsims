#' Calculate SSB0 at any given point in time 
#'
#' @param df Data frame with life history parameters 
#' @param time which timestep -idx - should ssb0 be calcualted (integer) 
#'
#' @return Returns ssbs0
#' @export
#'
#' @examples
#' 
calcSSB0 <- function(df, time = df$nyear){
  
  R0 <- exp(df$parms$logRinit)
  M <- df$M0[time]
  wage <- df$wage_ssb[,time]
  max_age <- df$age_maxage
  mat <- df$Matsel
  
  
  
  mage <- max_age # Max age
  age <- 0:(mage)
  nage <- length(age)

  if(length(M) == 1){
    M <- M*df$Msel
  }
  
    
  Mtmp <- rep(NA, nage)
  
  Mtmp[1:nage] <- cumsum(M)
  
  N0tmp <- rep(NA,nage)
  
  N0tmp[1] <- R0
  N0tmp[2:(nage-1)] = R0*exp(-Mtmp[1:(nage-2)])
  N0tmp[nage] =  R0*exp(-Mtmp[nage-1])/(1-exp(-Mtmp[nage]))
  
  N0 <- matrix(NA,nage)
  N0[1:(nage-1)] <- N0tmp[1:(nage-1)]
  N0[nage] <- N0tmp[nage]
  
  SSB0 <- sum(N0*wage*mat)

  
  
  return(SSB0)
}