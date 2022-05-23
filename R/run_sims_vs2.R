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
SDR <- 0.4
rho <- c(.01,.1, .5, .7, .9) # Try five different auto correlation coefficients 
#rho <- rep(.5, 5)
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


F0 <- 0.1
nruns <- 50
years <- 50


nruns <- 100
rho <- c(.9)#, 0.3, 0.5, .9)
set.seed(123)

df.in <- fn_sims(  models = c('linear'),
                   recLambda = c('noBOFF','BOFF'),
                   nruns = nruns,
                   years = years,
                   tau = tau,
                   Linf = Linf,
                   maxage = maxage,
                   K = K,
                   t0 = 0,
                   SDR = 0.2,
                   F0 = F0,
                   M = M0,
                   tau_sel = tau_sel,
                   egg.scale = egg.scale,
                   R0 = R0,
                   rho = rho,
                   Fpast = F0,
                   lambda.in = c(0.1, 0.5, .8),
                   recruitment = 'BH_R',
                   lambda.slope = .2,
                   mortality = 'constant',
                   fishing.type = 'AR',
                   recruitment.type = 'AR'
                   
)


df.in2 <- fn_sims(  models = c('linear'),
                   recLambda = c('noBOFF','BOFF'),
                   nruns = nruns,
                   years = years,
                   tau = tau,
                   Linf = Linf,
                   maxage = maxage,
                   K = K,
                   t0 = 0,
                   SDR = 1,
                   F0 = F0,
                   M = M0,
                   tau_sel = tau_sel,
                   egg.scale = egg.scale,
                   R0 = R0,
                   rho = rho,
                   Fpast = F0,
                   lambda.in = c(0.1, 0.5, .8),
                   recruitment = 'BH_R',
                   lambda.slope = .2,
                   mortality = 'constant',
                   fishing.type = 'AR',
                   recruitment.type = 'AR'
                   
)

ggplot(df.out[df.out$run == 1,], aes(x = propOld, y= rec, color = factor(lambda)))+geom_point()+geom_smooth()



#   if(i == 1){
df.out <- df.in[[1]]
df.N <- df.in[[3]]
df.save <- df.in[[2]]


df.out$lambda[is.na(df.out$lambda)] <- 0
df.save$lambda[is.na(df.out$lambda)] <- 0
df.N$lambda[is.na(df.out$lambda)] <- 0

#   if(i == 1){
df.out2 <- df.in2[[1]]
df.N2 <- df.in2[[3]]
df.save2 <- df.in2[[2]]


df.out2$lambda[is.na(df.out$lambda)] <- 0
df.save2$lambda[is.na(df.out$lambda)] <- 0
df.N2$lambda[is.na(df.out$lambda)] <- 0


df.plot <- rbind(df.out, df.out2)
df.plot$lambda[df.plot$lambda == 0] <- 1
df.plot$R0base <- 1000


# Plot some stuff 
prop.plot <- df.plot[df.plot$years>round(max(df.plot$years)/2),] %>% 
  group_by(model, run, Linf, SDR, K, tau, rho, lambda) %>% 
  dplyr::summarise( SSBcorRR0  = cor(SSBprop, recR0),
                    weightcorRR0  = cor(mweight, recR0),                                            
                    agecorRR0  = cor(mage, recR0),                       
                    SSBcorRtot  = cor(SSBprop, rtot),                         
                    weightcorRtot  = cor(mweight, rtot),                       
                    agecorRtot  = cor(mage, rtot),                      
                    SSBcorRdev  = cor(SSBprop, Rresid),                       
                    #weightcorRdev  = cor(mweight, Rresid),                     
                    agecorRdev  = cor(mage, Rresid),
                    SSBcorR = cor(SSBprop, rec),
                    magecorR = cor(mage, rec)
  ) %>% 
  pivot_longer(9:17)


dodge <- position_dodge(width = 0.8)


ggplot(prop.plot[prop.plot$name == 'SSBcorRdev',], aes(x = factor(SDR), y = value, color = factor(rho)))+geom_boxplot(width = 0.1,position = dodge)+
  geom_hline(aes(yintercept = 0), linetype = 2)+
  facet_grid(model ~ lambda)+theme_classic()


ggplot(prop.plot[prop.plot$name == 'SSBcorRtot',], aes(x = factor(SDR), y = value, color = factor(rho)))+geom_boxplot(width = 0.1,position = dodge)+
  geom_hline(aes(yintercept = 0), linetype = 2)+
  facet_grid(model ~ lambda)+theme_classic()



# Plot some other things 
df.save$SDR <- .2
df.save2$SDR <- 1

df.savetot <- rbind(df.save, df.save2) %>%  
  group_by(years, rho, M, lambda, model , SDR) %>% 
  summarise(Bio = median(SSB),
            Biomin = quantile(SSB, probs = 0.05),
            Biomax = quantile(SSB, probs = 0.95),
            Rec = median(R),
            Rmin = quantile(R, probs = 0.05),
            Rmax = quantile(R, probs = 0.95),
            R0 = median(R0_boff)
            )

df.savetot$lambda[is.na(df.savetot$lambda)] <- 1

p1 <- ggplot(df.savetot ,aes( x =years, y = Bio, color = factor(lambda) ))+geom_line()+facet_grid(SDR ~ rho)+theme_bw()
p2 <- ggplot(df.savetot ,aes( x =years, y = Rec, color = factor(lambda) ))+geom_line()+facet_grid(SDR ~ rho)+theme_bw()

p3 <- ggplot(df.savetot ,aes( x =years, y = R0, color = factor(lambda) ))+
  geom_line()+facet_grid(SDR ~ rho)+theme_bw()

ggsave(filename = 'figures/9-5/biotime.png', p1, width = 16, height = 12)
ggsave(filename = 'figures/9-5/rectime.png', p2, width = 16, height = 12)
ggsave(filename = 'figures/9-5/R0time.png', p3, width = 16, height = 12)



dodge <- position_dodge(width = 0.8)

p4 <- ggplot(prop.plot, 
       aes(x = factor(SDR), y = value, fill = factor(rho)))+geom_violin(position = dodge)+
  geom_boxplot(width = .2, position = dodge)+
  geom_hline(aes(yintercept = 0), linetype = 2)+
  facet_grid(lambda ~ name)+theme_classic()+coord_cartesian(ylim = c(-.4,.4))+
  scale_y_continuous('correlation')

ggsave(filename = 'figures/9-5/correlations.png', p4, width = 16, height = 12)


# Try something else 




