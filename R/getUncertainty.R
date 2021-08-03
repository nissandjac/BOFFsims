getUncertainty <- function(name,
                           data,
                           sdrep,
                           rep.values,
                           convers = NA){
  
  df <- data.frame(value = sdrep[rep.values == name,1])
  df$SE <- sdrep[rep.values == name,2]
  df$min <- df[,1]-2*df$SE
  df$max <- df[,1]+2*df$SE
  df$metric <- name
  
  
  
  if(dim(df)[1] == data$nyear){
    df$year <- 1:data$nyear
  }
  
  if(is.na(convers)){
    # Do nothing
  }else{
    df[1:4] <- exp(df[1:4])
  }
  
  return(df)  
}

