# Run MSE # 
library(reshape2)
library(gtools)
library(patchwork)
library(tidyverse)
library(TMB)

source('R/calcSSB0.R')
source('R/run_agebased_model_true_Catch.R')
source('R/load_data_seasons.R')


# Load the data frame with eggs 

eggs <- read.csv('data/fecundityEggSizeFemaleSize.csv')
eggs$weight <- 0.01*(eggs$FemaleSize_mm/10)^3 # Fix this later
# All eggs as a function of size 
ggplot(eggs, aes(x = weight, y = exp(lnFecundity), color = Species))+
  geom_point()+theme_bw()+theme(legend.position = 'none')+
  geom_smooth(method = 'lm')+
  scale_x_log10()+
  scale_y_log10()

# 
unique(eggs$Species)

# Just take the cod eggs most of the other fish are not really big fisheries 
cod <- eggs[eggs$Species == 'Gadus morhua',]

cod$weight <- 0.01*(cod$FemaleSize_mm/10)^3 # Fix the parameters for weight lenght to whatever here 



nruns <- 100
saving <- TRUE
seeds <- round(runif(n = nruns, min = 1,  max = 1e6))


df <- load_data_seasons(nseason = 4,
                        nyear = 100,# Set up parameters 
                        Linf = 150, 
                        maxage = 10,
                        K = 1, 
                        t0 = 0, SDR = .7,
                        fishing.type = 'constant',
                        mortality = 'constant',
                        alpha = 1e6,
                        beta = 2) # Specify parameters


lmegg <- lm(exp(lnFecundity) ~ weight, data = cod)

# Assume a 0 intercept 
intercept <- 0.
lmegg <- lm(I(exp(lnFecundity) - intercept) ~ 0 + weight, cod)
lmegglog <- lm(I(lnFecundity - intercept) ~ 0 + weight, cod)

# Try custom function 




wmodel <- seq(1,max(cod$weight), length.out = 100)

beta <- 2
eggs.hyp <- wmodel^2*exp(rnorm(length(wmodel), sd = 0.2))
eggs.lin <- (wmodel*beta*1e4)*exp(rnorm(length(wmodel), sd = 0.2))

plot(wmodel, eggs.hyp)
points(wmodel,eggs.lin, col = 'red')


plot((wmodel/eggs.hyp)[2:100])
plot((wmodel/eggs.lin)[2:100])


df.t <- data.frame(weight = rep(wmodel, 2), 
                   fecundity = c(
                     wmodel*lmegg$coefficients[1],
                     exp(wmodel*lmegglog$coefficients[1])
                   ),
                   model = rep(c('linear','log'), each = length(wmodel))
)
ggplot(df.t, aes(x = weight, y = fecundity, color = model))+geom_line()+theme_bw()


ggplot(cod, aes(x = weight, y = Fecundity_nOfEggs_per_female))+geom_point()+
  theme_classic()+
  geom_line(data = df.t, aes(x = weight,  y = fecundity, color = model))+
  scale_y_log10()


df.t$fecundity[df.t$fecundity < 0] <- NA

ggplot(cod, aes(x = weight, y = Fecundity_nOfEggs_per_female/weight))+geom_point()+
  theme_classic()+
  geom_line(data = df.t, aes(x = weight,  y = fecundity/weight, color = model))



tmprun <- run.agebased.true.catch(df, seed = seeds[1])


df$egg.size <- 



for(i in 1:nruns){
  set.seed(seeds[i])
  
  df <- load_data_seasons(nseason = 4,
                          nyear = 100,# Set up parameters 
                          Linf = 30, 
                          maxage = 10,
                          K = 1, 
                          t0 = 0, SDR = .7,
                          fishing.type = 'constant',
                          mortality = 'constant',
                          alpha = 1e6,
                          beta = 2) # Specify parameters 
  
  
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

