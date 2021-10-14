est_eggs <- function(x,y){
  
  
  
  xfun <- function(data,par){
    
    alpha <- par[[1]]
    
    if(length(par)>1){
      beta <- par[[2]]
      }else{
      beta <- 1
    }
    
    
    ml <- alpha*data[,1]^beta
    
    nll <- -sum(dnorm(ml, data[,2], sd = .1, log = TRUE))
    return(nll)
  }
  
  data <- cbind(x,y)
    
  
  xfun(data, par = list(100))
  
  
  
  par <- list(alpha = 0.1)
  
  result.lin <- optim(par = par, fn = xfun,
                  lower = 0.00001, upper = 1e6,
                  method = 'Brent',
                  data = data)
  
  
  ml.lin <- result.lin$par*x
  
  plot(x,ml.lin)
  points(x,y)
  
  print(paste('slope in linear model is = ',round(result.lin$par[1],digits = 3)))
  
  
  par <- list(alpha = 0.1,
                beta = 1)
    
    result <- optim(par = par, fn = xfun, data = data, 
                    method = 'L-BFGS-B',
                    lower = c(0.001,1),
                    upper = c(Inf,4))
    
    ml.hyper <- result$par[1]*x^result$par[2]
    
    print(paste('factor in hyper model is = ',round(result$par[1],digits = 3),
                ', exponent is calculated as =', round(result$par[2], digits = 3)))
   

  df <- data.frame(weight = rep(x, 3),
                   estimate = c(ml.lin, ml.hyper,y),
                   model = rep(c('linear','hyper','data'), each = length(x)))
    
    
  p1 <- ggplot(df[df$model != 'data',], aes(x = weight, y= estimate, color = model))+
    geom_line(size = 1.2)+
    geom_point(data = df[df$model == 'data',], alpha = 0.2)+
    theme_classic()+scale_y_continuous('fecundity',expand = c(0.01,0.01))+
    scale_x_continuous(expand = c(0.01,0.01))+
    coord_cartesian(xlim = c(0, max(df$weight)))
    
  print(p1)
  
  ls.parms <- list(alpha.lin = result.lin$par[[1]],
                   beta.lin = 1,
                   alpha.hyper = result$par[1],
                   beta.hyper = result$par[[2]])
  
  ls.return <- list(df = df,
                    parameters = ls.parms)
    
  return(ls.return)
 
  
  
   }
  
  
  
  
  
  
