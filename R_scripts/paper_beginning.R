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
source('R/getR0.R')

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


cod.df <- codest$df

p1 <- ggplot(cod.df %>% filter(model != 'data'), aes(x = weight, y = estimate, color = model))+geom_line()+
  geom_point(data=cod.df %>% filter(model == 'data'), alpha = .5)+scale_y_continuous('Fecundity')+
  scale_x_continuous('weight (g)')+theme_classic()
p1


# Do it for sardine as well

sard <- eggs %>% filter(Species =="Engraulis mordax")

sard$weight <- 0.01*(sard$FemaleSize_mm/10)^3 # Fix the parameters for weight lenght to whatever here 

sardest <- est_eggs(x = sard$weight,
                    y = sard$Fecundity_nOfEggs_per_female)


sard.df <- sardest$df

p2 <- ggplot(sard.df %>% filter(model != 'data'), aes(x = weight, y = estimate, color = model))+geom_line()+
  geom_point(data=sard.df %>% filter(model == 'data'), alpha = .5)+scale_y_continuous('Fecundity')+
  scale_x_continuous('weight (g)')+theme_classic()
p2


# Plot them together
pp.df <- rbind(cod.df %>% mutate(species = 'cod'), sard.df %>% mutate(species = 'sardine'))


p2 <- ggplot(pp.df %>% filter(model != 'data'), aes(x = weight, y = estimate, color = model))+geom_line()+
  geom_point(data=pp.df %>% filter(model == 'data'), alpha = .5)+scale_y_continuous('Fecundity')+
  facet_wrap(~species, scales = 'free', nrow = 1)+
  scale_x_continuous('weight (g)')+theme_classic()+theme(legend.position = c(0.2,0.8), legend.title = element_blank())
p2

ggsave(p2, file = 'egg_production_hyper_iso.png', width = 16, height = 8, units = 'cm')



Linf <- c(30, 60, 150)
K <- c(.7, 0.4, .2)
SDR <- c(1, 0.5, .4)
rho <- c(.1, .5, .8)
maxage <- c(5, 10, 20)
tau <- c(1.5,4,5)
M0 <- c(0.5, .4, .2)
tau_sel <- c(.5, 3, 3)
h <- c(0.7, .5, .4)
#egg.scale <- c(5e4, 3e5, 8e5)
R0 = c(10000,10000,10000)#*codest$parameters[['alpha.lin']]
recruitment = 'BH_steep'

nspecies <- 3 
nruns <- 10
seeds <- round(runif(nruns, 1, 1e6 ))
# Calculate Fmsy as a function of LH parameters 


for(n in 1:nruns){
  if(n == 1){
    start_time <- Sys.time()
  }
  for(i in 1:nspecies){
    
    Fmsy <- seq(0,M0[i]*2, length.out = 50)
    
    
    for(j in 1:length(Fmsy)){
      
      df <- load_data_seasons(nseason = 1,
                              nyear = 50,# Set up parameters 
                              Linf = Linf[i], 
                              maxage = maxage[i],
                              tau = tau[i],
                              K = K[i], 
                              h = h[i],
                              t0 = 0, 
                              tau_sel = tau_sel[i],
                              M= M0[i],
                              SDR = SDR[i], # Recruitment deviations - set to zero to calculate lambda
                              fishing.type = 'constant',
                              mortality = 'constant',
                              recruitment = recruitment,
                              negg = codest$parameters[['alpha.lin']],
                              eggbeta = codest$parameters[['beta.lin']],
                              F0 = Fmsy[j], # Set to zero to calc lambda
                              R0 = R0[i], 
                              seed = seeds[n]) # Specify parameters
      
      if(n == 1){
        R0[i] <- getR0(init = 10000, df)
        df$parms$logRinit <- log(R0[i])
      }
      
      
      
      tmp <- run.agebased.true.catch(df)
      
      
      if(i == 1 & j == 1 & n == 1){
        df.iso <- data.frame(yr = df$years, species = Linf[i], SSB = as.numeric(tmp$SSB), 
                             F0 = Fmsy[j], 
                             C = tmp$Catch, R = tmp$R.save, n = n)
      }else{
        df.iso <- rbind(df.iso,
                        data.frame(yr = df$years, species = Linf[i], 
                                   SSB = as.numeric(tmp$SSB), 
                                   F0 = Fmsy[j], 
                                   C = tmp$Catch, R = tmp$R.save, n = n))
      }
      
      
      
    }
  }
  
  for(i in 1:nspecies){
    
    Fmsy <- seq(0,M0[i]*2, length.out = 50)
    
    
    for(j in 1:length(Fmsy)){
      
      df <- load_data_seasons(nseason = 1,
                              nyear = 50,# Set up parameters 
                              Linf = Linf[i], 
                              maxage = maxage[i],
                              tau = tau[i],
                              K = K[i], 
                              h = h[i],
                              t0 = 0, 
                              tau_sel = tau_sel[i],
                              M= M0[i],
                              SDR = SDR[i], # Recruitment deviations - set to zero to calculate lambda
                              fishing.type = 'constant',
                              mortality = 'constant',
                              recruitment = recruitment,
                              negg = codest$parameters[['alpha.hyper']],
                              eggbeta = codest$parameters[['beta.hyper']],
                              F0 = Fmsy[j], # Set to zero to calc lambda
                              R0 = R0[i], 
                              seed = seeds[n]) # Specify parameters
      if(n == 1){
        R0[i] <- getR0(init = 10000, df)
        df$parms$logRinit <- log(R0[i])
      }
      
      
      tmp <- run.agebased.true.catch(df)
      
      
      if(i == 1 & j == 1 & n == 1){
        df.hyp <- data.frame(yr = df$years, species = Linf[i], 
                             SSB = as.numeric(tmp$SSB), 
                             F0 = Fmsy[j], C = tmp$Catch, R = tmp$R.save, n = n)
      }else{
        df.hyp <- rbind(df.hyp,
                        data.frame(yr = df$years, species = Linf[i], 
                                   SSB = as.numeric(tmp$SSB), 
                                   F0 = Fmsy[j], C = tmp$Catch,R = tmp$R.save, n = n))
      }
      
      
      
    }
  }
  if(n == nruns){
    end_time <- Sys.time()
    print(end_time - start_time)
  }
}  


df.iso <- df.iso %>% mutate(model = 'isometric')
df.hyp <- df.hyp %>% mutate(model = 'hyperallometric')


df.out <- rbind(df.hyp, df.iso)

df.out <- df.out %>% filter(yr %in% 40:50) %>% group_by(species, F0, model) %>% summarise(SSB = median(SSB),
                                                           SSBmin = quantile(SSB, probs = 0.05),
                                                           SSBmax = quantile(SSB, probs = 0.95),
                                                           Cm = median(C),
                                                           Cmin = quantile(C, probs = 0.05),
                                                           Cmax = quantile(C, probs = 0.95),
                                                           Rm = median(R),
                                                           Rmin = quantile(R, probs = 0.05),
                                                           Rmax = quantile(R, probs = 0.95)) %>% 
  group_by(species, model) %>% 
  mutate(scaleint = max(Cmax)) %>% 
  ungroup() %>% 
  mutate(Cm.scaled = Cm/scaleint,
         Cmin.scaled = Cmin/scaleint,
         Cmax.scaled = Cmax/scaleint) %>% select(-scaleint)

df.out$size <- 'small'
df.out$size[df.out$species == Linf[2]] <- 'medium'
df.out$size[df.out$species == Linf[3]] <- 'large'
df.out$size <- factor(df.out$size, levels = c('small','medium','large'))


MSY <- df.out %>% group_by(model, species, size) %>% summarise(MSY = max(Cm),
                                                               Fmsy = F0[which.max(Cm)])


p.2 <- ggplot(df.out, aes(x = F0, y = Cm.scaled, color = model, fill = model))+
  geom_line()+geom_ribbon(aes(ymin = Cmin.scaled, ymax = Cmax.scaled), linetype = 0, alpha = 0.2, show.legend = FALSE)+
 # geom_vline(data = MSY, aes(xintercept = Fmsy, color = model))+
  #geom_hline(data = MSY, aes(yintercept = MSY, color = model))+
  theme_classic()+
  facet_wrap(~size, scales = 'free_y', nrow = 3)+
  theme(legend.position = 'top', legend.title = element_blank())+
  scale_y_continuous('Catch')
p.2

ggsave(p.2, file = 'MSY_hyper_iso.png', width = 8, height = 12, units = 'cm')



## Plot MSY ## 
msy.p <- rbind(df.hyp, df.iso)
