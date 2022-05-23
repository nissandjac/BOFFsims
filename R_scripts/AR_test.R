# Test AR model # 
F0 <- 0.1
nruns <- 30
years <- 100


nruns <- 100
rho <- c(.1,.2 ,.5,.9)#, 0.3, 0.5, .9)

t0 = 0
SDR = 0.6
M <- .2
lambda.cut <- .9
lambda <- NA
R0 <- 1000
               
recruitment = 'BH_steep'
mortality = 'constant'
fishing.type = 'AR'
recruitment.type = 'AR'

cod <- eggs[eggs$Species == 'Gadus morhua',]

cod$weight <- 0.01*(cod$FemaleSize_mm/10)^3 # Fix the parameters for weight lenght to whatever here 

codest <- est_eggs(x = cod$weight,
                   y = cod$Fecundity_nOfEggs_per_female)

### Calculate Fmsy 


negg <- codest$parameters[['alpha.lin']]/egg.scale
eggbeta = codest$parameters[['beta.lin']]


  for(i in 1:length(rho)){         
    set.seed(123)
    
    df <- load_data_seasons(nseason = 1,
                        nyear = years,# Set up parameters 
                        Linf = Linf, 
                        maxage = maxage,
                        K = K, 
                        t0 = t0, 
                        M = M,
                        tau_sel = tau_sel,
                        tau = tau,
                        SDR = SDR, # Recruitment deviations 
                        fishing.type = fishing.type,
                        mortality = mortality,
                        alpha = alpha,
                        beta = beta,
                        recruitment = recruitment,
                        recruitment.type = recruitment.type,
                        negg = negg,
                        eggbeta = eggbeta,
                        F0 = F0[1],
                        rhoR = rho[i],
                        R0 = R0,
                        lambda = lambda,
                        lambda.cut = lambda.cut) # Specify parameters



tmprun <- run.agebased.true.catch(df)



if(i == 1){
  df.save <- data.frame(R = tmprun$R.save,
                        time = 1:df$tEnd,
                        rho = rho[i])
}else{
  df.save <- rbind(df.save, data.frame(R = tmprun$R.save,
                                       time = 1:df$tEnd,
                                       rho = rho[i]))
}

}


ggplot(df.save, aes(x = time, y = R, color = factor(rho)))+geom_point()+geom_line()+theme_classic()


