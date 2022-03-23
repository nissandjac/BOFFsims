### Run bulk sims ###
require('scam')
require('tidyverse')
require('patchwork')
source('R/fn_sims.R')
source('R/calcSSB0.R')
source('R/run_agebased_model_true_Catch.R')
source('R/load_data_seasons.R')
source('R/est_eggs.R')
source('R/plotRecruitment.R')
source('R/getEquilibrium.R')
source('R/runScenarios.R')

eggs <- read.csv('data/fecundityEggSizeFemaleSize.csv')
eggs$weight <- 0.01*(eggs$FemaleSize_mm/10)^3 # Fix this later
# All eggs as a function of size 

eggs <- eggs %>% group_by(Species) %>% 
  mutate(relweight = weight/max(weight),                                     
         releggs = Fecundity_nOfEggs_per_female/max(Fecundity_nOfEggs_per_female)
  )


x <- eggs[is.na(eggs$relweight) == 0,]$relweight
y <- eggs[is.na(eggs$relweight) == 0,]$releggs


# all relative eggs 
#parms <- est_eggs(x,y)



# Just take the cod eggs most of the other fish are not really big fisheries 
cod <- eggs[eggs$Species == 'Gadus morhua',]

cod$weight <- 0.01*(cod$FemaleSize_mm/10)^3 # Fix the parameters for weight lenght to whatever here 

codest <- est_eggs(x = cod$weight,
                   y = cod$Fecundity_nOfEggs_per_female)


# Life history 

Linf <- 50
K <-  0.4
SDR <- 0.5
rho <- c(.01, .5, .99) # Try three different auto correlation coefficients 
maxage <- 10
tau <- 3
M0 <- .4
tau_sel <- 2
egg.scale <- 1e5
R0 = 1000
recruitment = 'BH_R'

nspecies <- 1
F0 <- seq(0.1,1, length.out = 5)


# Calculate Fmsy as a function of LH parameters 
Fmsy <- seq(0,2, length.out = 50)

Catch <- data.frame(
  F0 = rep(Fmsy, nspecies),
  species = rep(Linf, each = length(Fmsy)*nspecies),
  Catch = NA
)


for(i in 1:nspecies){
  for(j in 1:length(Fmsy)){
    
    df <- load_data_seasons(nseason = 1,
                            nyear = 100,# Set up parameters 
                            Linf = Linf[i], 
                            maxage = maxage[i],
                            tau = tau[i],
                            K = K[i], 
                            t0 = 0, 
                            tau_sel = tau_sel[i],
                            M= M0[i],
                            SDR = 0, # Recruitment deviations - set to zero to calculate lambda
                            fishing.type = 'constant',
                            mortality = 'constant',
                            recruitment = recruitment,
                            negg = codest$parameters[['alpha.lin']]/egg.scale[i],
                            eggbeta = codest$parameters[['beta.lin']],
                            F0 = Fmsy[j], # Set to zero to calc lambda
                            R0 = R0) # Specify parameters
    
    
    tmp <- run.agebased.true.catch(df)
    
    
    Catch[Catch$species == Linf[i] & Catch$F0 == Fmsy[j],]$Catch <- tmp$Catch[df$tEnd]
    
    
  }
}


ggplot(Catch, aes(x = F0, y = Catch, color = as.factor(species)))+geom_line()+theme_classic()+
  facet_wrap(~species, scales = 'free_y')





for(i in 1:nspecies){
  
  
  df.in <- fn_sims(tau = tau[i],
                   Linf = Linf[i],
                   maxage = maxage[i],
                   K = K[i],
                   t0 = 0,
                   SDR = SDR[i],
                   F0 = F0,
                   M = M0[i],
                   tau_sel = tau_sel[i],
                   egg.scale = egg.scale[i],
                   R0 = 1000,
                   rho = rho[i],
                   recruitment = 'BH_R',
                   lambda.slope = .7,
                   mortality = 'constant',
                   fishing.type = 'constant',
                   recruitment.type = 'AR'
                   
  )
  
  if(i == 1){
    df.out <- df.in
  }else{
    df.out <- rbind(df.out, df.in)
  }
  
  
}


# remove the ones with zero recruitment 
df.plot <- df.out[df.out$rec>1,]



p1 <- ggplot(df.plot[df.plot$years > 50 ,], aes(x = mweight, y = SSBprop, color = model))+geom_point()+
  facet_wrap(~Linf, scales = 'free')+theme_bw()


p2 <- ggplot(df.plot[df.plot$years > 50 ,], aes(x = SSBprop, y = Rdev, color = model))+geom_point()+
  facet_wrap(~Linf, scales = 'free')+theme_bw()+geom_smooth(method = 'lm')


p3 <- ggplot(df.plot[df.plot$years > 50 ,], aes(x = SSBprop, y = recR0, color = model))+geom_point()+
  facet_wrap(~Linf, scales = 'free')+theme_bw()+geom_smooth(method = 'lm')

p1/p2/p3

prop.plot <- df.plot %>% group_by(model, run, Linf, SDR, K, tau) %>% summarise(SSBcorRR0  = cor(SSBprop, recR0),
                                                                               weightcorRR0  = cor(mweight, recR0),
                                                                               agecorRR0  = cor(mage, recR0),
                                                                               SSBcorRtot  = cor(SSBprop, rtot),
                                                                               weightcorRtot  = cor(mweight, rtot),
                                                                               agecorRtot  = cor(mage, rtot),
                                                                               SSBcorRdev  = cor(SSBprop, Rresid),
                                                                               weightcorRdev  = cor(mweight, Rresid),
                                                                               agecorRdev  = cor(mage, Rresid)) %>% 
  pivot_longer(7:15)


lims <- max(abs(prop.plot$value))


p1 <-  ggplot(prop.plot, aes(x = model, y = value, fill = model))+geom_violin()+geom_boxplot(width = 0.1)+
  facet_grid(factor(Linf)~name)+scale_x_discrete()+
  theme_bw()+theme(axis.text.x = element_blank())+coord_cartesian(ylim = c(-lims,lims))+geom_hline(aes(yintercept = 0), linetype = 2)


ggsave('figures/all_species_violin.png', p1, width = 16, height = 16)
