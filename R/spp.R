

K <- 1000
t <- 50
r <- 0.2 


B <- rep(0, t)


F0 <- seq(0,1, length.out = 30)
C <- B <- matrix(NA, t, length(F0))
B[1,] <- K

for(k in 1:length(F0)){
  for(i in 1:(t-1)){
  
    B[i+1,k] <- B[i,k]+r*B[i,k]*(1-B[i,k]/K)-B[i,k]*F0[k] 
    C[i,k] <- B[i,k]*F0[k]
  }

  
}

plot(F0, C[49,], type = 'l')
max(C[49,], na.rm = TRUE)
