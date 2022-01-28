# Run egg trial thingie # 
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
source('R/runScenarios.R')
# Load the data frame with eggs 

eggs <- read.csv('data/fecundityEggSizeFemaleSize.csv')
eggs$weight <- 0.01*(eggs$FemaleSize_mm/10)^3 # Fix this later
# All eggs as a function of size 

eggs <- eggs %>% group_by(Species) %>% mutate(relweight = weight/max(weight),
                                              releggs = Fecundity_nOfEggs_per_female/max(Fecundity_nOfEggs_per_female)
)


x <- eggs[is.na(eggs$relweight) == 0,]$relweight
y <- eggs[is.na(eggs$relweight) == 0,]$releggs


# all relative eggs 
parms <- est_eggs(x,y)



# Just take the cod eggs most of the other fish are not really big fisheries 
cod <- eggs[eggs$Species == 'Gadus morhua',]

cod$weight <- 0.01*(cod$FemaleSize_mm/10)^3 # Fix the parameters for weight lenght to whatever here 

codest <- est_eggs(x = cod$weight,
                   y = cod$Fecundity_nOfEggs_per_female)

png('hyper_normal_eggs.png', width = 18, height = 18, units = 'cm', res = 400)
codest <- est_eggs(x = cod$weight,
                   y = cod$Fecundity_nOfEggs_per_female)
dev.off()


ggplot(codest$df[codest$df$model != 'data',], 
       aes(x =  weight, y = estimate/weight, color = model))+geom_line()+theme_bw()

nruns <- 100
saving <- TRUE
seeds <- round(runif(n = nruns, min = 1,  max = 1e6))



# Punch in life history parameters (these are standard in run scenario)
tau <- 5
Linf <- 150
maxage <- 10
K <- 0.6
t0 <- 0
SDR <- 0.5
F0 <- 0.2
R0 <- 1000
M <- 0.4
recruitment <- 'BH_R'
lambda.slope <- .6



ls.plot <- runScenarios(models = c('linear','hyper'),
                        recLambda = c('noBOFF','BOFF'),
                        runLambda = FALSE,
                        lambda.in = .3, # has to be a fraction
                        egg.df = codest, 
                        lambda.slope = lambda.slope,
                        SDR = SDR,
                        F0 = F0)

# Run all the models 

p1 <- ggplot(ls.plot[[1]][ls.plot[[1]]$years > 50,], aes(x = years, y = S, color = model, fill = model))+geom_line()+
  geom_ribbon(aes(ymin = Smin, ymax = Smax), show.legend = FALSE, alpha = 0.1, linetype = 0)+theme_classic()

p2 <- ggplot(ls.plot[[1]][ls.plot[[1]]$years > 50,], aes(x = years, y = C, color = model, fill = model))+geom_line()+
  geom_ribbon(aes(ymin = Cmin, ymax = Cmax), show.legend = FALSE, alpha = 0.1, linetype = 0)+theme_classic()

p3 <- ggplot(ls.plot[[1]][ls.plot[[1]]$years > 50,], aes(x = years, y = Rec, color = model, fill = model))+geom_line()+
  geom_ribbon(aes(ymin = Recmin, ymax = Rmax), show.legend = FALSE, alpha = 0.1, linetype = 0)+theme_classic()


# 

# p4 <- ggplot(ls.plot[[1]][ls.plot[[1]]$years > 50,], aes(x = years, y = Rtot, color = model, fill = model))+
#   geom_line()+theme_classic()
# p4


png(file = 'stuff.png', res = 400, height = 16, width = 16, units = 'cm')
p1/p2/p3
dev.off()

# Run them over fishing mortalities 

Frun <- seq(0,0.7, length.out = 20)
lambda <- seq(0.2, 0.9, length.out = 20)


for(i in 1:length(Frun)){
#  for(j in 1:length(lambda)){  
  
    ls.plot <- runScenarios(models = c('linear','hyper'),
                          recLambda = c('noBOFF','BOFF'),
                          runLambda = FALSE,
                          lambda.in = .4, # has to be a fraction
                          egg.df = codest, 
                          lambda.slope = lambda.slope,
                          SDR = SDR,
                          F0 = Frun[i])
  
  
  if(i == 1){
    df.save <- ls.plot[[2]]
    df.save$F0 <- Frun[i]
  #  df.save$lambda <- lambda[j]
  }else{
    
    tmp <- ls.plot[[2]]
    tmp$F0 <- Frun[i]
   # tmp$lambda <- lambda[j]
    
    df.save <- rbind(df.save, tmp)
    
  }
  
#  } 
}

# Summarise the df by fishing mortality 

df.sum <- df.save %>% 
  group_by(years, model, F0) %>% 
  summarise(SSB = median(SSB),
            C = median(Catch),
            R = median(R)) %>% pivot_longer(c('SSB','C','R'))


p5 <- ggplot(df.sum[df.sum$years == 100,], aes(x = F0,  y = value, color = model))+geom_line()+
  theme_classic()+facet_wrap(~name, scales = 'free')

p5

# 
# png('withF0.png', )

#
