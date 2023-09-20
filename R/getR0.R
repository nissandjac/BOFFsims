getR0 <- function(init = 10000, 
                  df){
  
  nage <- df$nage
  N0 <- rep(0, nage)
  

  Mage <- cumsum(df$M0)
  
  
  N0[1] <- init
  N0[2:(nage-1)] = init*exp(-Mage[1:(nage-2)])
  N0[nage] =  init*exp(-Mage[nage-1])/(1-exp(-df$M0[nage]))
  
  
  R_out <- sum(N0 * df$Matsel * df$egg.size)
  
 return(R_out) 
}