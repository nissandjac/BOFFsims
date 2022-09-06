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
SDR <- 0.
rho <- c(.01,.1, .5, .7, .9) # Try five different auto correlation coefficients 
#rho <- rep(.5, 5)
maxage <- 10
tau <- 3
M0 <- .4
tau_sel <- 2
egg.scale <- 1e5
R0 = 100

recruitment = 'BH_steep'
recruitment.type <- 'AR'
fishing.type <- 'constant'
nspecies <- 1
mortality <- 'constant'
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

Fmsy = Catch$F0[which.max(Catch$Catch)]

years <- 100
F0 <- c(seq(0.1,Fmsy*1.5, length.out = floor(years/2)),seq(Fmsy*1.5,Fmsy/2, length.out = ceiling(years/2)))
t0 <- 0
models = c('linear')
recLambda = c('noBOFF','BOFF')
recruitment <- 'BH_steep'
mortality = 'constant'
fishing.type = 'constant'
recruitment.type = 'AR'

nruns <- 300
rho <- c(0, 0.3, 0.5, .9)
set.seed(123)
lambda.in <- c(1.5)
lambda.cut <- .9
M <- 0.2


set.seed(21353)

ls.plot <- runScenarios(models = models,
                        recLambda = recLambda,
                        nruns = nruns, 
                        years = years,
                        Fpast = Fpast,
                        runLambda = FALSE,
                        lambda.in = lambda.in,
                        rho = rho,
                        egg.df = codest, 
                        egg.scale = egg.scale,
                        lambda.cut = lambda.cut,
                        SDR = .4,
                        h = .2,
                        SDF = 0,
                        F0 = F0,
                        maxage = maxage,
                        K = K, 
                        Linf = Linf,
                        t0 = t0,
                        tau = tau,
                        tau_sel = tau_sel,
                        M = M,
                        mortality = mortality,
                        fishing.type = fishing.type,
                        recruitment.type = recruitment.type
)


# 
pfrac <- ggplot(ls.plot$df.save[ls.plot$df.save$years > 20,], aes(x = F0,y = xfrac,
                                                         color = model))+geom_point(alpha =.05)+geom_smooth(se = FALSE)+
  scale_y_continuous('Fraction of old spawners (x)')+theme_classic()+theme(legend.position = c(.8,.8),
                                                                             legend.background=element_blank())


p1 <- ggplot(ls.plot$df.save[ls.plot$df.save$years > 20,], aes(x = xfrac,y = R,color = model))+
  geom_point(alpha =.05)+geom_smooth(size = 1.5)+
  theme_classic()+scale_y_continuous('Recruitment')+scale_x_continuous('Fraction of \nold spawners (x)')+
  theme(legend.position = 'none')


plot1 <- p1+pfrac+ plot_annotation(tag_levels = 'a')

plot1


ggsave('figboff_F0.png',plot1 , width = 16, height = 10, units = 'cm')
# Calculate fraction young and old based on mikaels calcs 

df.N <- ls.plot$df.N[ls.plot$df.N$age > 0,]
df.sum <- ls.plot$df.sum

df.N$old <- 'young'
df.N$old[df.N$age >= (df$tau+2)] <- 'old'


xx <- df.N %>% group_by(years, F0, rho, lambda, model, old, run ) %>% summarise(SSB = sum(N*weight*mat),
                                                                                Rdev = mean(Rdev))

# Exclude recruits from all calculations 
xx.tot <- df.N %>% group_by(years, rho, lambda, model, run ) %>% 
  summarise(SSB = sum(N*weight*mat),
            Rdev = mean(Rdev),
            mage = weighted.mean(age, N*mat),
            mweight = weighted.mean(weight, N*mat))

xxOld <- xx[xx$old == 'old',]
xxYoung <- xx[xx$old == 'young',]

# Sort the 
fracOld <- data.frame(oldFrac = xxOld$SSB/(xxOld$SSB+xxYoung$SSB),
                      F0 = xxOld$F0,
                      years = xxOld$years,
                      lambda = xxOld$lambda,
                      rho = xxOld$rho,
                      Rdev = xxOld$Rdev,
                      mSSBage  = xx.tot$mage,
                      mSSBweight = xx.tot$mweight,
                      model = xxOld$model,
                      run = xxOld$run)


## Make a correlation data frame for Mikael ## 

df.cor <- fracOld[fracOld$years > 20,]  %>% group_by(F0, lambda, rho, run, model) %>% 
  mutate(Rdev_frac = cor(oldFrac, Rdev),
         Rdev_mage = cor(oldFrac, mSSBage),
         Rdev_mweight = cor(oldFrac, mSSBweight)
         ) #%>%  select(years, F0, lambda, rho, Rdev_frac, mage_frac, mweight_frac, run)
    
fracOld$F0 <- round(fracOld$F0, digits = 2)


p1 <- ggplot(fracOld[fracOld$run == 1 & fracOld$years > 20,], aes(x = years, y = oldFrac, group = model, color = factor(F0)))+
  geom_line()+facet_wrap(~F0)+
  geom_hline(aes(yintercept = mean(fracOld$oldFrac)), linetype = 2)+
  theme(legend.position = 'none')+theme_classic()+theme(legend.position = 'none')+
  scale_y_continuous('fraction of old spawners')


p2 <- ggplot(fracOld[fracOld$run == 1 & fracOld$years > 20,], aes(x = F0, y = oldFrac, group = model, color = factor(F0)))+
  geom_point()+facet_wrap(~model)+
  geom_hline(aes(yintercept = mean(fracOld$oldFrac)), linetype = 2)+
  theme(legend.position = 'none')+theme_classic()+theme(legend.position = 'none')+
  scale_y_continuous('')


p1+p2

ggsave('figures/F0_impact_onspawner.png', p1+p2, width = 16, height = 12, units = 'cm')


# Plot the correlation plots # 
df.cor.plot <- df.cor %>% pivot_longer(c(Rdev_frac, Rdev_mage, Rdev_mweight)) 

dodge <- position_dodge(.5)

ggplot(df.cor.plot, aes(x = name, y = value, color = model))+geom_violin(position = dodge)+
  geom_boxplot(width = .2,position = dodge)+theme_bw()


p3 <- ggplot(df.cor, aes(x = oldFrac, y = Rdev))+geom_point()+geom_smooth()+scale_x_continuous('')+theme_classic()+
  scale_y_continuous('recruitment\ndeviations')+theme_classic()
p4 <- ggplot(df.cor, aes(x = oldFrac, y = mSSBage))+geom_point()+geom_smooth()+
  scale_x_continuous('Fraction of\nold spawners')+theme_classic()+
  scale_y_continuous('mean age')
p5 <- ggplot(df.cor, aes(x = oldFrac, y = mSSBweight))+geom_point()+geom_smooth()+theme_classic()+
  scale_y_continuous('mean weight')+scale_x_continuous('')


p3 + p4 + p5

ggsave('figures/cors.png', p3+p4+p5, width = 16, height = 12, units = 'cm')

p5 <- ggplot(df.cor, aes(x = oldFrac, y = Rdev))+geom_point()+geom_smooth()
p6 <- ggplot(df.cor, aes(x = mSSBage, y = Rdev))+geom_point()+geom_smooth()
p7 <- ggplot(df.cor, aes(x = mSSBweight, y = Rdev))+geom_point()+geom_smooth()

p5 + p6 + p7

ggsave('figures/rdev_calcs.png', p3+p4+p5, width = 16, height = 12, units = 'cm')


# Now calculate the post recruitment SR relationship 

R.df <- ls.plot$df.save

# Change this to whatever model combinations are running 
modelruns <- paste(R.df$model, R.df$rho, sep = '-')
# 

R.df$modelruns <- modelruns
R.df$residuals <- NA
R.df$Rpred <- NA

for(i in 1:nruns){
  for(j in 1:length(unique(modelruns))){
    
    
    dftmp <- R.df[R.df$run == i & R.df$modelruns == unique(modelruns)[j],]  
    
    Rtmp <- scam(log(R) ~ s(SSB, k = 20, bs = 'mpd', m = 2) +  
                   offset(log(SSB)), 
                 family=gaussian(link="identity"), data = dftmp, optimizer="nlm",sp=0.01)
    
    
    # Calculate the residuals (anomalies) 
    dftmp$Rpred <- exp(predict(Rtmp, newdata = data.frame(SSB = dftmp$SSB)))
    dftmp$residuals <- log(dftmp$Rpred)-log(dftmp$R)
    
    
    # Add the residuals to the R.df data frame 
    R.df[R.df$run == i & R.df$modelruns == unique(modelruns)[j],]$residuals <- dftmp$residuals
    R.df[R.df$run == i & R.df$modelruns == unique(modelruns)[j],]$Rpred <- dftmp$Rpred
    
  }
}

# Pick 10 random  runs for better plots 
idx <- round(runif(n = 10, min = 1, max = 100))

pR <- ggplot(R.df[R.df$run %in% idx,], aes(x = SSB/1e6, y = R, color = factor(rho)))+
  geom_point(alpha =0.1)+
  geom_line(aes(y = Rpred))+
  facet_grid(model~run)+theme_classic()+theme(legend.position = 'top')+
  scale_x_continuous('SSB (1000 kg)',breaks = c(0, 1, 2))
pR

ggsave('recruitment_scam.png',pR, width = 16, height = 10, units = 'cm')


# Redo the correlation plot # 
df.N <- ls.plot$df.N[ls.plot$df.N$age > 0,]
df.N$Rdev_scam <- rep(R.df$residuals, each = max(df$age))
df.sum <- ls.plot$df.sum

df.N$old <- 'young'
df.N$old[df.N$age >= (df$tau+2)] <- 'old'

xx <- df.N %>% group_by(years, rho, lambda, model, old, run ) %>% 
  summarise(SSB = sum(N*weight*mat), 
            Rdev = mean(Rdev),
            Rdev_scam = mean(Rdev_scam))

# # Exclude recruits from all calculations 
xx.tot <- df.N %>% group_by(years, rho, lambda, model, run ) %>%
  summarise(SSB = sum(N*weight*mat),
            Rdev = mean(Rdev),
            mage = weighted.mean(age, N*mat),
            mweight = weighted.mean(weight, N*mat))

xxOld <- xx[xx$old == 'old',]
xxYoung <- xx[xx$old == 'young',]

# Sort the 
fracOld <- data.frame(oldFrac = xxOld$SSB/(xxOld$SSB+xxYoung$SSB),
                      #F0 = xxOld$F0,
                      years = xxOld$years,
                      lambda = xxOld$lambda,
                      rho = xxOld$rho,
                      Rdev = xxOld$Rdev,
                      Rdev_scam = xxOld$Rdev_scam,
                      mSSBage  = xx.tot$mage,
                      mSSBweight = xx.tot$mweight,
                      model = xxOld$model,
                      run = xxOld$run)


# Violin plots # 

p_dev1 <- ggplot(fracOld[fracOld$years>15,], aes(x = oldFrac, y = Rdev_scam, color = model))+geom_point(alpha = .1)+
  geom_smooth(method = 'lm', size = 1.3)+theme_classic()+
  scale_x_continuous('fraction of \nold spawners')+
  scale_y_continuous('Rec dev from scam')+
  theme(legend.position = 'none')

p_dev2 <- ggplot(fracOld[fracOld$years>15,], aes(x = oldFrac, y = Rdev, color = model))+geom_point(alpha = .1)+
  geom_smooth(method = 'lm', size = 1.3, se = FALSE)+theme_classic()+
  scale_x_continuous('fraction of \nold spawners')+
  scale_y_continuous('True rec devs')+
  theme(legend.position = c(.8,.9))

p_dev1+p_dev2


# Correlation coefficients to match figure 4 

df.cor <- fracOld[fracOld$years > 15,]  %>% group_by(lambda, rho, run, model) %>% 
  mutate(Rdev_old = cor(oldFrac, Rdev),
         Rdev_old_scam = cor(oldFrac, Rdev_scam),
         Rdev_mage = cor(mSSBage,Rdev),
         Rdev_mage_scam = cor(mSSBage,Rdev_scam),
         Rdev_mweight = cor(mSSBweight, Rdev),
         Rdev_mweight_scam = cor(mSSBweight, Rdev_scam),
  ) %>% pivot_longer(c(Rdev_old,Rdev_old_scam, Rdev_mage,Rdev_mage_scam,
                       Rdev_mweight,Rdev_mweight_scam), values_to = 'correlation')

dodge <- position_dodge(0.2)


df.cor$model[df.cor$model == 'linear-BOFF'] <- 'BOFF'
df.cor$model[df.cor$model == 'linear-noBOFF'] <- 'Standard'

pcor <- ggplot(df.cor, aes(x = name, y = correlation, color = model))+geom_violin(position = dodge)+
  geom_boxplot(width = .05,position = dodge)+
  scale_y_continuous('correlation coefficient with deviations')+
  facet_wrap(~rho, ncol = 1)+geom_hline(aes(yintercept = 0), linetype = 2, color = 'black')+
  theme_classic()+scale_x_discrete('', labels = c('%old spawners\n(true)',
                                                  '%old spawners\n(scam)',
                                                  'mean age \n(true)','mean age \n(scam)', 
                                                  'mean weight\n(true)','mean weight \n(scam)'))+
  theme(axis.text.x = element_text(angle = 45, vjust = 0, hjust=0),
        legend.position = 'top',
        legend.title = element_blank())

pcor


#windows(width = 10/cm(1), height = 16/cm(1))

ggsave('correlation_plot.png',pcor, width = 10, height = 16, units = 'cm')


