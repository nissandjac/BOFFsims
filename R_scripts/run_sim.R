# Run MSE from the RAM stocks # 
library(reshape2)
library(gtools)
library(patchwork)
library(tidyverse)
library(TMB)

source('R/calcSSB0.R')
source('R/run_agebased_model_true_Catch.R')
source('R/load_data_seasons.R')
source('R/load_data_future.R')

nyear <- 100
nruns <- 100
set.seed(1234)
saving <- TRUE
seeds <- round(runif(n = nruns, min = 1,  max = 1e6))

df <- load_data_seasons(nseason = 4,
                        nyear = nyear,# Set up parameters 
                        Linf = 30, 
                        maxage = 10,
                        K = 1, 
                        t0 = 0, SDR = .7,
                        fishing.type = 'AR',
                        mortality = 'AR') # Specify parameters 


df.save <- data.frame(years = rep(df$years, nruns),
                      SSB = NA,
                      R = NA,
                      Catch = NA, 
                      run = rep(1:nruns, each = length(df$years)),
                      model = 'standard')

for(i in 1:nruns){
  set.seed(seeds[i])
  
  df <- load_data_seasons(nseason = 4,
                          nyear = nyear,# Set up parameters 
                          Linf = 30, 
                          maxage = 10,
                          K = 1, 
                          t0 = 0, SDR = .7,
                          fishing.type = 'AR',
                          mortality = 'AR') # Specify parameters 
  
  tmprun <- run.agebased.true.catch(df, seed = seeds[i])
  
  df.save[df.save$run == i,]$SSB <- tmprun$SSB
  df.save[df.save$run == i,]$R <- tmprun$R
  df.save[df.save$run == i,]$Catch <- tmprun$Catch
  
}


df.sum <- df.save %>% group_by(years, model) %>% 
  summarise(S = mean(SSB),
            Rec = mean(R),
            C = mean(Catch),
            )



df.save.boff <- data.frame(years = rep(df$years, nruns),
                      SSB = NA,
                      R = NA,
                      Catch = NA, 
                      run = rep(1:nruns, each = length(df$years)),
                      model = 'boff')

for(i in 1:nruns){
  set.seed(seeds[i])
  
  df <- load_data_seasons(nseason = 4,
                          nyear = nyear,# Set up parameters 
                          Linf = 30, 
                          maxage = 10,
                          K = 1, 
                          t0 = 0, SDR = .7,
                          fishing.type = 'AR',
                          mortality = 'AR',
                          beta = 1.3) # Specify parameters 
  
  tmprun <- run.agebased.true.catch(df, seed = seeds[i])
  
  df.save.boff[df.save.boff$run == i,]$SSB <- tmprun$SSB
  df.save.boff[df.save.boff$run == i,]$R <- tmprun$R
  df.save.boff[df.save.boff$run == i,]$Catch <- tmprun$Catch
  
}


df.sum.boff <- df.save.boff %>% group_by(years, model) %>% 
  summarise(S = mean(SSB),
            Rec = mean(R),
            C = mean(Catch),
  )




# Plot boff and not boff
df.plot <- rbind(df.save, df.save.boff)
df.sumplot <- rbind(df.sum, df.sum.boff)

ggplot(df.plot, aes(x= years, y = SSB, group = as.factor(run), color = model))+geom_line(alpha = 0.1, size = .5)+theme_bw()+theme(legend.position = 'none')+
  geom_line(data = df.sumplot, aes(y = S, group = NA), size = 2)



ggplot(df.sumplot, aes(x = years, y = S, color = model))+geom_line()+theme_bw()+
  geom_hline(aes(yintercept = tmprun$SSB_0))

