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
lambda.slope <- .7


df <- load_data_seasons(nseason = 1,
                        nyear = 100,# Set up parameters 
                        Linf = Linf, 
                        maxage = maxage,
                        tau = tau,
                        K = K, 
                        t0 = t0, 
                        M= M,
                        SDR = 0, # Recruitment deviations - set to zero to calculate lambda
                        fishing.type = 'constant',
                        mortality = 'constant',
                        recruitment = recruitment,
                        negg = codest$parameters[['alpha.lin']],
                        eggbeta = codest$parameters[['beta.lin']],
                        F0 = 0, # Set to zero to calc lambda
                        R0 = R0) # Specify parameters




Fin <- seq(0,1, length.out = 15)

ls.plot <- runScenarios(models = c('linear','hyper'),
                        recLambda = c('noBOFF','BOFF'),
                        runLambda = FALSE,
                        lambda.in = .4,
                        egg.df = codest, 
                        lambda.slope = lambda.slope,
                        SDR = SDR,
                        F0 = Fin)

# Run all the models 

p1 <- ggplot(ls.plot[[1]][ls.plot[[1]]$years > 50 & ls.plot[[1]]$F0 == Fin[4],], aes(x = years, y = S, color = model, fill = model))+geom_line()+
  geom_ribbon(aes(ymin = Smin, ymax = Smax), show.legend = FALSE, alpha = 0.1, linetype = 0)+theme_classic()+theme(legend.position = 'top')

p2 <- ggplot(ls.plot[[1]][ls.plot[[1]]$years > 50 & ls.plot[[1]]$F0 == Fin[4],], aes(x = years, y = C, color = model, fill = model))+geom_line()+
  geom_ribbon(aes(ymin = Cmin, ymax = Cmax), show.legend = FALSE, alpha = 0.1, linetype = 0)+theme_classic()+theme(legend.position = 'none')

p3 <- ggplot(ls.plot[[1]][ls.plot[[1]]$years > 50 & ls.plot[[1]]$F0 == Fin[4],], aes(x = years, y = Rec, color = model, fill = model))+geom_line()+
  geom_ribbon(aes(ymin = Recmin, ymax = Rmax), show.legend = FALSE, alpha = 0.1, linetype = 0)+theme_classic()+theme(legend.position = 'none')

 p4 <- ggplot(ls.plot[[1]][ls.plot[[1]]$years > 50 & ls.plot[[1]]$F0 == Fin[4],], aes(x = years, y = Rtot, color = model, fill = model))+
   geom_line()+theme_classic()+theme(legend.position = 'none')
# p4


p1/p2/p3/p4


png('figures/boff_stuff_combined.png', res = 400, width = 20, height = 15, units = 'cm')
p1/p2/p3/p4
dev.off()

p5 <- ggplot(ls.plot[[1]][ls.plot[[1]]$years == 99,], aes(x = F0,y = Rec, color = model, fill = model))+geom_line(size = 1.1)+
  geom_ribbon(aes(ymin = Recmin, ymax = Rmax), show.legend = FALSE, alpha = 0.1, linetype = 0)+theme_classic()+theme(legend.position = 'top')
p6 <- ggplot(ls.plot[[1]][ls.plot[[1]]$years == 99,], aes(x = F0,y = S, color = model, fill = model))+geom_line(size = 1.1)+
  geom_ribbon(aes(ymin = Smin, ymax = Smax), show.legend = FALSE, alpha = 0.1, linetype = 0)+theme_classic()+theme(legend.position = 'none')
p7 <- ggplot(ls.plot[[1]][ls.plot[[1]]$years == 99,], aes(x = F0,y = C, color = model, fill = model))+geom_line()+
  geom_ribbon(aes(ymin = Cmin, ymax = Cmax), show.legend = FALSE, alpha = 0.1, linetype = 0)+theme_classic()+theme(legend.position = 'none')
p8 <- ggplot(ls.plot[[1]][ls.plot[[1]]$years == 99,], aes(x = F0,y = Rtot, color = model))+geom_line()+theme_classic()+theme(legend.position = 'none')

png('figures/boff_stuff_F0.png', res = 400, width = 20, height = 15, units = 'cm')
p5/p6/p7/p8
dev.off()



# Make Mikaels figures # 

df.N <- ls.plot[[3]]

df.N$SSB <- df.N$N*df.N$weight*df.N$mat


## Summarise the data frame to  plot it 
df.N$old <- NA
pold <- tau

df.N$old[df.N$age < pold] <- 'young'
df.N$old[df.N$age >= pold] <- 'old'



df.Nsum <- df.N %>% 
  group_by(years, F0, model, old, run) %>% 
  summarise(N = sum(N),
            Catch = sum(Catch),
            SSB=  sum(SSB),
            Rdev = mean(Rdev)) # The Rdev mean is not a mean (it's the same for all ) 

# Get median weight weighted by numbers

df.wSum <- df.N %>% 
  group_by(years, F0, model, run) %>%
  summarise(mWeight = weighted.mean(weight, N),
            mAge = weighted.mean(age, N))
  



df.propOld <- data.frame(propOld = df.Nsum$N[df.Nsum$old == 'young']/df.Nsum$N[df.Nsum$old == 'old'],
                         SSBprop = df.Nsum$SSB[df.Nsum$old == 'young']/df.Nsum$SSB[df.Nsum$old == 'old'],
                         years = df.Nsum$years[df.Nsum$old == 'old'],
                         model = df.Nsum$model[df.Nsum$old == 'old'],
                         run = df.Nsum$run[df.Nsum$old == 'old'],
                         meanWeight = df.wSum$mWeight,
                         meanAge = df.wSum$mAge,
                         R = ls.plot[[2]]$R,
                         RR0 = ls.plot[[2]]$R/exp(df$parms$logRinit),
                         Rtot = ls.plot[[2]]$Rtot
                         )



median(df.propOld$propOld)
median(df.propOld$SSBprop) # approx 50% with age 6

# Summarise the data 


df.propSum <- df.propOld %>% 
  group_by( years, model) %>% 
 summarise(propN = median(propOld),
           propSSB = median(SSBprop),
           mWeight = median(meanWeight),
           mAge = median(meanAge)) %>% 
  pivot_longer(cols = 3:6)



ggplot(df.propSum[df.propSum$years > 50,], aes(x = years, y= value, color = model))+geom_line()+facet_wrap(~name, scales = 'free_y')

# Try the correlation coefficients 


propplot <- df.propOld %>% group_by(model, run) %>% summarise(SSBcor  = cor(SSBprop, RR0),
                                                              weightcor  = cor(meanWeight, RR0),
                                                              agecor  = cor(meanAge, RR0)) %>% pivot_longer(3:5)


p1 <- ggplot(propplot, aes(x = model, y = value, group = model, fill = model))+geom_violin()+facet_wrap(~name)+theme_bw()+geom_hline(aes(yintercept = 0))+
  geom_boxplot(width = 0.2)+theme(legend.position = 'top',
                                  axis.text.x=element_blank())+scale_x_discrete('')+
  scale_y_continuous('correlation coefficient (SSB ~ R/R0)')



propplot <- df.propOld %>% group_by(model, run) %>% summarise(SSBcor  = cor(SSBprop, Rtot),
                                                         weightcor  = cor(meanWeight, Rtot),
                                                         agecor  = cor(meanAge, Rtot)) %>% pivot_longer(3:5)

p2 <- ggplot(propplot, aes(x = model, y = value, group = model, fill = model))+
  geom_violin()+
  facet_wrap(~name)+
  theme_bw()+
  geom_hline(aes(yintercept = 0))+
  geom_boxplot(width = 0.2)+
  theme(legend.position = 'none',axis.text.x = element_text(angle = 50, vjust = .5))+
  scale_x_discrete('')+
  scale_y_continuous('correlation coefficient (SSB ~ #eggs)')



p1/p2







