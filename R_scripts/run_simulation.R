# Run MSE # 
library(reshape2)
library(gtools)
library(patchwork)
library(tidyverse)
library(TMB)

source('R/calcSSB0.R')
source('R/run_agebased_model_true_Catch.R')
source('R/load_data_seasons.R')
source('R/est_eggs.R')
source('R/plotRecruitment.R')
source('R/getEquilibrium.R')

# Load the data frame with eggs 

eggs <- read.csv('data/fecundityEggSizeFemaleSize.csv')
eggs$weight <- 0.01*(eggs$FemaleSize_mm/10)^3 # Fix this later
# All eggs as a function of size 

eggs <- eggs %>% group_by(Species) %>% mutate(relweight = weight/max(weight),
                                              releggs = Fecundity_nOfEggs_per_female/max(Fecundity_nOfEggs_per_female)
)

cod <- eggs[eggs$Species == 'Gadus morhua',]

cod$weight <- 0.01*(cod$FemaleSize_mm/10)^3 # Fix the parameters for weight lenght to whatever here 

codest <- est_eggs(x = cod$weight,
                   y = cod$Fecundity_nOfEggs_per_female)


nruns <- 100
seeds <- round(runif(n = nruns, min = 1,  max = 1e6))


# Calculate max number of eggs assuming 1 recruit distribution 

Linf <- 150
maxage <- 12
K <- 0.4
t0 <- 0
SDR <- 0.5
F0 <- 0.5
R0 <- 1000
recruitment <- 'BH_R'


# Calculate linear rep as a function of F

nruns <- 100 


for(i in 1:nruns){

  set.seed(seeds[i])
  
  df <- load_data_seasons(nseason = 1,
                          nyear = 100,# Set up parameters 
                          Linf = Linf, 
                          maxage = maxage,
                          K = K, 
                          t0 = t0, 
                          SDR = SDR, # Recruitment deviations 
                          fishing.type = 'constant',
                          mortality = 'constant',
                          alpha = alpha,
                          beta = beta,
                          recruitment = recruitment,
                          recruitment.type = 'AR',
                          negg = codest$parameters[['alpha.lin']],
                          eggbeta = codest$parameters[['beta.lin']],
                          F0 = F0,
                          R0 = R0) # Specify parameters
  
  tmprun <- run.agebased.true.catch(df)
  
  # Iterate over fishing mortalities
  ls.lin <- Iteratesims_F(df, Fin = round(seq(0,.5, length.out = 20), digits = 2)) # Run over 30 fishing mortalities
  
  if(i == 1){
    df.out <- ls.lin[[1]]
    df.out$run <- paste('run', i, sep = '-')
  }else{
    tmp <- ls.lin[[1]]
    tmp$run <- paste('run', i, sep = '-')
    
    df.out <- rbind(df.out,tmp)
  }

}


for(i in 1:nruns){
  
  set.seed(seeds[i])
  
  df <- load_data_seasons(nseason = 1,
                          nyear = 100,# Set up parameters 
                          Linf = Linf, 
                          maxage = maxage,
                          K = K, 
                          t0 = t0, 
                          SDR = SDR, # Recruitment deviations 
                          fishing.type = 'constant',
                          mortality = 'constant',
                          alpha = alpha,
                          beta = beta,
                          recruitment = recruitment,
                          recruitment.type = 'AR',
                          negg = codest$parameters[['alpha.hyper']],
                          eggbeta = codest$parameters[['beta.hyper']],
                          F0 = F0,
                          R0 = R0) # Specify parameters
  
  tmprun <- run.agebased.true.catch(df)
  
  # Iterate over fishing mortalities
  ls.lin <- Iteratesims_F(df, Fin = round(seq(0,.5, length.out = 20), digits = 2),
                          model = 'hyper') # Run over 30 fishing mortalities
  
  if(i == 1){
    df.out.hyper <- ls.lin[[1]]
    df.out.hyper$run <- paste('run', i, sep = '-')
  }else{
    tmp <- ls.lin[[1]]
    tmp$run <- paste('run', i, sep = '-')
    
    df.out.hyper <- rbind(df.out.hyper,tmp)
  }
  
}


df.tot <- rbind(df.out,df.out.hyper)

## Summarize 


df.sum <- df.tot %>% group_by(years, model,F0) %>% 
  summarise(S = median(SSB),
            Rec = median(R),
            C = median(Catch),
            Rtot = median(Rtot),
            Smin = quantile(SSB, probs = 0.05),
            Smax = quantile(SSB, probs = 0.95),
            Recmin = quantile(R, probs = 0.05),
            Rmax = quantile(R, probs = 0.95),
            Cmin = quantile(Catch, probs = 0.05),
            Cmax = quantile(SSB, probs = 0.95)
  )

# 
# ggplot(df.sum[df.sum$years > 30,], aes(x = years, y = Rtot, color = model))+geom_line()+
#   facet_wrap(~F0, scales = 'free_y'))
#   theme_bw()

ggplot(df.sum[df.sum$years > 30,], aes(x = years, y = S, color = model))+geom_line()+
  facet_wrap(~F0, scales = 'free_y')+theme_bw()+
  geom_ribbon(aes(ymin = Smin, ymax = Smax, fill = model), alpha = 0.1, linetype = 0)

ggplot(df.sum[df.sum$years > 30,], aes(x = years, y = C, color = model))+geom_line()+
  facet_wrap(~F0, scales = 'free_y')+theme_bw()+
  geom_ribbon(aes(ymin = Cmin, ymax = Cmax, fill = model), alpha = 0.1, linetype = 0)



# Run simulation with SDR 

nruns <- 100
SDRin <- round(seq(0,1, length.out = 20), digits = 2)


for(i in 1:nruns){
  
  set.seed(seeds[i])
  
  df <- load_data_seasons(nseason = 1,
                          nyear = 100,# Set up parameters 
                          Linf = Linf, 
                          maxage = maxage,
                          K = K, 
                          t0 = t0, 
                          SDR = SDR, # Recruitment deviations 
                          fishing.type = 'constant',
                          mortality = 'constant',
                          alpha = alpha,
                          beta = beta,
                          recruitment = recruitment,
                          recruitment.type = 'AR',
                          negg = codest$parameters[['alpha.lin']],
                          eggbeta = codest$parameters[['beta.lin']],
                          F0 = F0,
                          R0 = R0) # Specify parameters
  
  tmprun <- run.agebased.true.catch(df)
  
  # Iterate over fishing mortalities
  ls.lin <- Iteratesims_SDR(df, SDRin = SDRin) # Run over 30 fishing mortalities
  
  if(i == 1){
    df.out <- ls.lin[[1]]
    df.out$run <- paste('run', i, sep = '-')
  }else{
    tmp <- ls.lin[[1]]
    tmp$run <- paste('run', i, sep = '-')
    
    df.out <- rbind(df.out,tmp)
  }
  
}


for(i in 1:nruns){
  
  set.seed(seeds[i])
  
  df <- load_data_seasons(nseason = 1,
                          nyear = 100,# Set up parameters 
                          Linf = Linf, 
                          maxage = maxage,
                          K = K, 
                          t0 = t0, 
                          SDR = SDR, # Recruitment deviations 
                          fishing.type = 'constant',
                          mortality = 'constant',
                          alpha = alpha,
                          beta = beta,
                          recruitment = recruitment,
                          recruitment.type = 'AR',
                          negg = codest$parameters[['alpha.hyper']],
                          eggbeta = codest$parameters[['beta.hyper']],
                          F0 = F0,
                          R0 = R0) # Specify parameters
  
  tmprun <- run.agebased.true.catch(df)
  
  # Iterate over fishing mortalities
  ls.lin <- Iteratesims_SDR(df,  SDRin =SDRin,
                          model = 'hyper') # Run over 30 fishing mortalities
  
  if(i == 1){
    df.out.hyper <- ls.lin[[1]]
    df.out.hyper$run <- paste('run', i, sep = '-')
  }else{
    tmp <- ls.lin[[1]]
    tmp$run <- paste('run', i, sep = '-')
    
    df.out.hyper <- rbind(df.out.hyper,tmp)
  }
  
}


df.tot <- rbind(df.out,df.out.hyper)

## Summarize 


df.sum <- df.tot %>% group_by(years, model,SDR) %>% 
  summarise(S = median(SSB),
            Rec = median(R),
            C = median(Catch),
            Rtot = median(Rtot),
            Smin = quantile(SSB, probs = 0.05),
            Smax = quantile(SSB, probs = 0.95),
            Recmin = quantile(R, probs = 0.05),
            Rmax = quantile(R, probs = 0.95),
            Cmin = quantile(Catch, probs = 0.05),
            Cmax = quantile(SSB, probs = 0.95)
  )

# 
# ggplot(df.sum[df.sum$years > 30,], aes(x = years, y = Rtot, color = model))+geom_line()+
#   facet_wrap(~F0, scales = 'free_y'))
#   theme_bw()

ggplot(df.sum[df.sum$years > 30,], aes(x = years, y = S, color = model))+geom_line()+
  facet_wrap(~SDR, scales = 'free_y')+theme_bw()+
  geom_ribbon(aes(ymin = Smin, ymax = Smax, fill = model), alpha = 0.1, linetype = 0)

ggplot(df.sum[df.sum$years > 30,], aes(x = years, y = C, color = model))+geom_line()+
  facet_wrap(~SDR, scales = 'free_y')+theme_bw()+
  geom_ribbon(aes(ymin = Cmin, ymax = Cmax, fill = model), alpha = 0.1, linetype = 0)






