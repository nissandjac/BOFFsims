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
nruns <- 30
years <- 100


nruns <- 100
rho <- c(.1,.2 ,.5,.9)#, 0.3, 0.5, .9)
set.seed(123)

df.in <- fn_sims(  nruns = nruns,
                   years = years,
                   tau = tau,
                   Linf = Linf,
                   maxage = maxage,
                   K = K,
                   t0 = 0,
                   SDR = 0.6,
                   F0 = F0,
                   M = M0,
                   tau_sel = tau_sel,
                   egg.scale = egg.scale,
                   R0 = R0,
                   rho = rho,
                   Fpast = F0,
                   recruitment = 'BH_R',
                   lambda.slope = .9,
                   mortality = 'constant',
                   fishing.type = 'AR',
                   recruitment.type = 'AR'
                   
  )
  
#   if(i == 1){
    df.out <- df.in[[1]]
    df.N <- df.in[[3]]
    df.save <- df.in[[2]]
  
    
    
ggplot(df.out[grep('linear', df.out$model),] , aes(x = mage, y = rec, color = model))+geom_point()+
  theme_classic()#+geom_smooth(method = 'lm', se = FALSE)
    

# Linear only 
ggplot(df.save, aes(x = years, y = R, color = model))+geom_point()+scale_y_log10()
ggplot(df.out, aes(x = years, y = mage, color = model))+geom_point()+scale_y_log10()+geom_smooth()
#   }else{
#     df.out <- rbind(df.out, df.in[[1]])
#     df.N <- rbind(df.N, df.in[[3]])
#     df.save <- rbind(df.save, df.in[[2]])
#   }
#   
#   
# }
# 

# remove the ones with zero recruitment 
df.plot <- df.out[df.out$rec>1,]



p1 <- ggplot(df.plot, aes(x = mweight, y = SSBprop, color = model))+geom_point(alpha = .2)+
  facet_wrap(~rho, scales = 'free')+theme_bw()+theme(legend.position = 'none')+geom_smooth(method = 'lm')
p1
# 
# 
# p2 <- ggplot(df.plot, aes(x = SSBprop, y = Rdev, color = model))+geom_point(alpha = .2)+
#   facet_wrap(~rho, scales = 'free')+theme_bw()+geom_smooth(method = 'lm')+theme(legend.position = 'none')
# 
# 
# p3 <- ggplot(df.plot, aes(x = SSBprop, y = recR0, color = model))+geom_point(alpha = .2)+
#   facet_wrap(~rho, scales = 'free')+theme_bw()+geom_smooth(method = 'lm')+theme(legend.position = 'none')
# 
# p1/p2/p3

# Look at the time series before correlating 
df.p <- df.save[df.save$model %in% c('linear-BOFF', 'linear-noBOFF'),]
df.p1 <- df.out[df.out$model %in% c('linear-BOFF', 'linear-noBOFF'),]


ggplot(df.p[df.save$run == 30,], aes(x = years, y = R, color = model))+geom_line()#+facet_wrap(~rho)


ggplot(df.p1, aes(x = mweight , y = rec, color = model))+geom_point()+geom_smooth()



# See if seeds are working correctly 

prop.plot <- df.plot[df.plot$years>round(max(df.plot$years)/2),] %>% group_by(model, run, Linf, SDR, K, tau, rho) %>% 
  dplyr::summarise( SSBcorRR0  = cor(SSBprop, recR0),
                    weightcorRR0  = cor(mweight, recR0),                                            
                    agecorRR0  = cor(mage, recR0),                       
                    SSBcorRtot  = cor(SSBprop, rtot),                         
                    weightcorRtot  = cor(mweight, rtot),                       
                    agecorRtot  = cor(mage, rtot),                      
                    SSBcorRdev  = cor(SSBprop, Rresid),                       
                    weightcorRdev  = cor(mweight, Rresid),                     
                    agecorRdev  = cor(mage, Rresid)
                    ) %>% 
  pivot_longer(8:16)

df.plot[df.plot$years>round(max(df.plot$years)/2),] %>% group_by(model, run, Linf, SDR, K, tau, rho) %>% 
  ggplot(aes(x = mage, y = Rresid))+geom_point(alpha =.2)+geom_smooth(method = 'lm')+theme_classic()


df.plot[df.plot$years>round(max(df.plot$years)/2),] %>% group_by(model, run, Linf, SDR, K, tau, rho) %>% 
  ggplot(aes(x = mweight, y = Rresid))+geom_point(alpha =.2)+geom_smooth(method = 'lm')+theme_classic()


df.plot[df.plot$years>round(max(df.plot$years)/2),] %>% group_by(model, run, Linf, SDR, K, tau, rho) %>% 
  ggplot(aes(x = mweight, y = Rresid))+geom_point(alpha =.2)+geom_smooth(method = 'lm')+theme_classic()


df.plot[df.plot$years>round(max(df.plot$years)/2),] %>% group_by(model, run, Linf, SDR, K, tau, rho) %>% 
  ggplot(aes(x = mweight, y = Rresid))+geom_point(alpha =.2)+geom_smooth(method = 'lm')+theme_classic()



lims <- max(abs(prop.plot$value))

p1 <-  ggplot(prop.plot, aes(x = model, y = value, fill = model))+geom_violin()+geom_boxplot(width = 0.1)+
  facet_grid(factor(rho)~name)+scale_x_discrete()+
  theme_bw()+theme(axis.text.x = element_blank())+coord_cartesian(ylim = c(-lims,lims))+geom_hline(aes(yintercept = 0), linetype = 2)


p2 <-  ggplot(prop.plot, aes(x = factor(rho), y = value, fill = factor(rho)))+geom_violin()+geom_boxplot(width = 0.1)+
  facet_grid(model~name)+scale_x_discrete()+
  theme_bw()+theme(axis.text.x = element_blank())+coord_cartesian(ylim = c(-lims,lims))+geom_hline(aes(yintercept = 0), linetype = 2)
p2



p1

ggsave('figures/rho_violin_nofishing.png', p1, width = 16, height = 16)
ggsave('figures/rhomodel_violin_nofishing.png', p2, width = 12, height = 16)

# Plot some stock recruitment 
idx <- c(1, 15, 25, 30, 50)
# Time series 
dfsave.tot <- df.save %>% 
  group_by(years, model,rho) %>% 
  dplyr::summarise(SSBmean = median(SSB), 
            SSBmin = quantile(SSB, 0.05),
            SSBmax = quantile(SSB, 0.95),
            Rmean = median(R), 
            Rmin = quantile(R, 0.05),
            Rmax = quantile(R, 0.95),
            Catchmean = median(Catch), 
            Catchmin = quantile(Catch, 0.05),
            Catchmax = quantile(Catch, 0.95))

p1 <- ggplot(dfsave.tot, aes(x = years, y = Rmean))+
  geom_line(data = df.save[df.save$run %in% idx,], aes(y = R, color = factor(run)), size = .6, alpha = .5)+
  geom_line(size = 1.2)+theme_classic()+facet_grid(rho~model)+
  geom_ribbon(aes(ymin = Rmin, ymax = Rmax), alpha = .2, linetype = 0)+
  theme(legend.position = 'none')
p1


p2 <- ggplot(dfsave.tot, aes(x = years, y = SSBmean))+
  geom_line(data = df.save[df.save$run %in% idx,], aes(y = SSB, color = factor(run)), size = .6, alpha = .5)+
  geom_line(size = 1.2)+theme_classic()+facet_grid(rho~model)+
  geom_ribbon(aes(ymin = SSBmin, ymax = SSBmax), alpha = .2, linetype = 0)+
  theme(legend.position = 'none')

p2

ggsave(p1, file  ='figures/recruitment_time.png', width = 20, height = 16, units = 'cm')
ggsave(p2, file = 'figures/SSB_time.png', width = 20, height = 16, units = 'cm')
# Use the scam fit for a range of models # 

models <- unique(df.save$model)
years <- unique(df.save$years)
nyears <- length(years)

df.save$residuals <- NA
df.save$Rscam <- NA
df.save$SSBx <- NA
df.save$SRtrue <- NA

# Scam calculation takes quite a while for many stocks 

for(i in 1:nruns){
  for(j in 1:length(models)){
    for(k in 1:length(rho)){
    
    dftmp <- df.save[df.save$run == i & df.save$model == models[j] & df.save$rho == rho[k],] 
    
    Ftmp <- scam(log(R+.001) ~ s(SSB, k = 20, bs = 'mpd', m = 2) +  
                   offset(log(SSB+.001)), 
                 family=gaussian(link="identity"), data = dftmp,optimizer="nlm",sp=0.01)
    
    
    # Calculate the residuals (anomalies) 
    dftmp$Rpred <- exp(predict(Ftmp, newdata = data.frame(SSB = dftmp$SSB)))
    dftmp$residuals <- log(dftmp$R)-log(dftmp$Rpred)
    
    
    # Add the residuals to the R.df data frame 
    df.save[df.save$run == i & df.save$model == models[j] & df.save$rho == rho[k],]$residuals <- dftmp$residuals
    df.save[df.save$run == i & df.save$model == models[j] & df.save$rho == rho[k],]$Rscam <- dftmp$Rpred
    
    Rtot <- dftmp$Rtot

    # Regulate R0 based on the lambda function 
    R0 <- dftmp$R0_boff
    
    #df.save[df.save$run == i & df.save$model == models[j] & df.save$rho == rho[k],]$SRtrue <- SR
    }
  }


}


p3 <- ggplot(df.save, aes(x = SSB,  y = Rscam, color = factor(run)))+geom_line()+geom_line()+facet_grid(rho ~ model)+
  theme_classic()+geom_point(aes(y = R), alpha = 0.2)+geom_line(aes(y = SR), col = 'black')+
  theme(legend.position = 'none')
p3

ggsave(p3, filename = 'figures/SR_relations.png', width = 20, height = 16, units = 'cm')


