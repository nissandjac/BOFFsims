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
SDR <- 0.2
rho <- c(.01,.1, .5, .7, .99) # Try five different auto correlation coefficients 
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
  facet_wrap(~species, scales = 'free_y')+theme(legend.position = 'none')


F0 <- 0.2


for(i in 1:length(rho)){
  
  
  df.in <- fn_sims(tau = tau,
                   Linf = Linf,
                   maxage = maxage,
                   K = K,
                   t0 = 0,
                   SDR = SDR,
                   F0 = F0,
                   M = M0,
                   tau_sel = tau_sel,
                   egg.scale = egg.scale,
                   R0 = R0,
                   rho = rho[i],
                   recruitment = 'BH_R',
                   lambda.slope = .7,
                   mortality = 'constant',
                   fishing.type = 'AR',
                   recruitment.type = 'AR'
                   
  )
  
  if(i == 1){
    df.out <- df.in[[1]]
    df.N <- df.in[[3]]
    df.save <- df.in[[2]]
  }else{
    df.out <- rbind(df.out, df.in[[1]])
    df.N <- rbind(df.N, df.in[[3]])
    df.save <- rbind(df.save, df.in[[2]])
  }
  
  
}


# remove the ones with zero recruitment 
df.plot <- df.out[df.out$rec>1,]



p1 <- ggplot(df.plot[df.plot$years > 50 ,], aes(x = mweight, y = SSBprop, color = model))+geom_point()+
  facet_wrap(~rho, scales = 'free')+theme_bw()


p2 <- ggplot(df.plot[df.plot$years > 50 ,], aes(x = SSBprop, y = Rdev, color = model))+geom_point()+
  facet_wrap(~rho, scales = 'free')+theme_bw()+geom_smooth(method = 'lm')


p3 <- ggplot(df.plot[df.plot$years > 50 ,], aes(x = SSBprop, y = recR0, color = model))+geom_point()+
  facet_wrap(~rho, scales = 'free')+theme_bw()+geom_smooth(method = 'lm')

p1/p2/p3

prop.plot <- df.plot %>% group_by(model, run, Linf, SDR, K, tau, rho) %>% summarise(SSBcorRR0  = cor(SSBprop, recR0),
                                                                               weightcorRR0  = cor(mweight, recR0),
                                                                               agecorRR0  = cor(mage, recR0),
                                                                               SSBcorRtot  = cor(SSBprop, rtot),
                                                                               weightcorRtot  = cor(mweight, rtot),
                                                                               agecorRtot  = cor(mage, rtot),
                                                                               SSBcorRdev  = cor(SSBprop, Rresid),
                                                                               weightcorRdev  = cor(mweight, Rresid),
                                                                               agecorRdev  = cor(mage, Rresid)) %>% 
  pivot_longer(8:16)


lims <- max(abs(prop.plot$value))


p1 <-  ggplot(prop.plot, aes(x = model, y = value, fill = model))+geom_violin()+geom_boxplot(width = 0.1)+
  facet_grid(factor(rho)~name)+scale_x_discrete()+
  theme_bw()+theme(axis.text.x = element_blank())+coord_cartesian(ylim = c(-lims,lims))+geom_hline(aes(yintercept = 0), linetype = 2)

p1

ggsave('figures/rho_violin.png', p1, width = 16, height = 16)


# Plot some stock recruitment 
idx <- c(1, 15, 25, 30, 50)
# Time series 

p2 <- ggplot(df.save[df.save$run %in% idx,], aes(x = years, y = R, color = factor(run)))+
  geom_line()+theme_classic()+facet_grid(rho~model)
p2


p3 <- ggplot(df.save[df.save$run %in% idx,], aes(x = years, y = SSB, color = factor(run)))+
  geom_line()+theme_classic()+facet_grid(rho~model)
p3


p2 <- ggplot(df.save[df.save$run %in% idx,], aes(x = years, y = Rtot, color = factor(run)))+
  geom_line()+theme_classic()+facet_grid(rho~model)
p2

p2 <- ggplot(df.save[df.save$run %in% idx,], aes(x = SSB, y = R, color = factor(run)))+
  geom_point()+theme_classic()+facet_grid(rho~model)
p2




# Use the scam fit for a range of models # 






















