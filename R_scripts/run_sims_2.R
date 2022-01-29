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
tau <- 2
Linf <- 30
maxage <- 5
K <- 0.8
t0 <- 0
SDR <- 1
F0 <- 0.2
R0 <- 1000
M <- 0.8
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


## Calculate Fmsy ###

Fmsy <- seq(.01, 3, length.out = 50)
Catch <- NA

for(i in 1:length(Fmsy)){
  
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
                          negg = codest$parameters[['alpha.lin']]/1e4,
                          eggbeta = codest$parameters[['beta.lin']],
                          F0 = Fmsy[i], # Set to zero to calc lambda
                          R0 = R0) # Specify parameters
  
  
 xx <- run.agebased.true.catch(df) 
 Catch[i] <- xx$Catch[df$tEnd]  
}

plot(Fmsy,Catch, type ='l')


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
            Rdev = mean(Rdev)) %>% arrange(run, old,model) # The Rdev mean is not a mean (it's the same for all ) 

# Get median weight weighted by numbers

df.wSum <- df.N %>% 
  group_by(years, F0, model ,run) %>%
  summarise(mWeight = weighted.mean(weight, N),
            mAge = weighted.mean(age, N)) %>% arrange(run, model)



# Rearrange ls.plot[22]

R.df <- ls.plot[[2]] %>% arrange(run, years,model)


df.propOld <- data.frame(propOld = df.Nsum$N[df.Nsum$old == 'old']/(df.Nsum$N[df.Nsum$old == 'young']+df.Nsum$N[df.Nsum$old == 'old']),
                         SSBprop = df.Nsum$SSB[df.Nsum$old == 'old']/(df.Nsum$SSB[df.Nsum$old == 'young']+df.Nsum$SSB[df.Nsum$old == 'old']),
                         years = df.Nsum$years[df.Nsum$old == 'old'],
                         model = df.Nsum$model[df.Nsum$old == 'old'],
                         run = df.Nsum$run[df.Nsum$old == 'old'],
                         Rdev = df.Nsum$Rdev[df.Nsum$old == 'old'],
                         mweight = df.wSum$mWeight,
                         mage = df.wSum$mAge,
                         rec = R.df$R,
                         M = R.df$M,
                         F0 = R.df$F0,
                         recR0 = R.df$R/exp(df$parms$logRinit),
                         rtot = R.df$Rtot
)

print(median(df.propOld$SSBprop))

# df.propSum <- df.propOld %>% 
#   group_by( years, model) %>% 
#  summarise(propN = median(propOld),
#            propSSB = median(SSBprop),
#            mWeight = median(meanWeight),
#            mAge = median(meanAge)) %>% 
#   pivot_longer(cols = 3:6)



# ggplot(df.propSum[df.propSum$years > 50,], aes(x = years, y= value, color = model))+geom_line()+facet_wrap(~name, scales = 'free_y')

# Try the correlation coefficients 

prop.plot <- df.propOld %>% group_by(model, run) %>% summarise(SSBcorRR0  = cor(SSBprop, recR0),
                                                               weightcorRR0  = cor(mweight, recR0),
                                                               agecorRR0  = cor(mage, recR0),
                                                               SSBcorRtot  = cor(SSBprop, rtot),
                                                               weightcorRtot  = cor(mweight, rtot),
                                                               agecorRtot  = cor(mage, rtot),
                                                               SSBcorRdev  = cor(SSBprop, Rdev),
                                                               weightcorRdev  = cor(mweight, Rdev),
                                                               agecorRdev  = cor(mage, Rdev)) %>%
  pivot_longer(3:11)


lims <- max(abs(prop.plot$value))


p1 <-   ggplot(prop.plot, aes(x = model, y = value, fill = model))+geom_violin()+geom_boxplot(width = 0.1)+
    facet_wrap(~name)+scale_x_discrete()+scale_y_continuous('correlation coefficient')+
    theme_bw()+theme(axis.text.x = element_blank())+coord_cartesian(ylim = c(-lims,lims))+geom_hline(aes(yintercept = 0), linetype = 2)

p1# 

png('figures/codlikeviolin.png', width = 16, height = 16, res = 400, units = 'cm')
p1
dev.off()






