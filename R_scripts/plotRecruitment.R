### Plot the four different recruitment functions ### 

  R0 <- 1000
  
  # BH first 
  alpha <- 1
  
  Rtot <- seq(1,1000, length.out = 150) 


  R_bh<- (alpha*Rtot)/(1+beta*Rtot)
  
  R_steep <- (4*h*R_0[space]*SSB[yr,space]/
          (SSB_0*(1-h)+ SSB[yr,space]*(5*h-1)))
  
  R_ricker <- alpha*Rtot *exp(-beta * Rtot)*exp(-0.5*df$b[yr]*SDR^2+Ry)#
  
  R_eggs <- R0*Rtot/(Rtot+R0)

  
  plot(Rtot,R_eggs, type = 'l')  

  df.plot -   