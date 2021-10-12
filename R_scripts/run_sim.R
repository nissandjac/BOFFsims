# Run MSE # 
library(reshape2)
library(gtools)
library(patchwork)
library(tidyverse)
library(TMB)

source('R/calcSSB0.R')
source('R/run_agebased_model_true_Catch.R')
source('R/load_data_seasons.R')
source('R/est_eggs.R')

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


# ggplot(eggs, aes(x = relweight, y = releggs, group = Species, color = Species))+
#   geom_point()+theme_bw()+theme(legend.position = 'none')+
#   geom_smooth(method = 'glm',se = FALSE)#+
# 
# 

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



alpha <- 0.05
beta <- 1e-4

SSBtest <- seq(1, 1/beta*1/alpha, length.out = 1000)
bhmodel <- plotRecruitment('Ricker', S = SSBtest, df = list(alpha = alpha, beta = beta))
plot(SSBtest,bhmodel, type = 'l')
lines(rep(1/beta, 100), y = seq(0,max(bhmodel), length.out =100))


df <- load_data_seasons(nseason = 1,
                        nyear = 200,# Set up parameters 
                        Linf = 150, 
                        maxage = 10,
                        h = 0.4,
                        K = .4, 
                        t0 = 0, 
                        SDR = 0, # Recruitment deviations 
                        fishing.type = 'constant',
                        mortality = 'constant',
                        alpha = alpha,
                        beta = beta,
                        recruitment = 'Ricker',
                        negg = codest$parameters[['alpha.lin']],
                        eggbeta = codest$parameters[['beta.lin']],
                        F0 = 0) # Specify parameters

df$egg.size <- df$egg.size/10

#df$parms$logRinit <- log(1)
tmp <- run.agebased.true.catch(df)

plot(tmp$SSB, type ='l', log = 'y')


plot(SSBtest,bhmodel, type ='l', xlim = c(0,max(tmp$Rtot.save)))
points(tmp$Rtot.save, tmp$R.save)



df.save <- data.frame(years = rep(df$years, nruns),
                           SSB = NA,
                           R = NA,
                           Rtot = NA,
                           Catch = NA, 
                           run = rep(1:nruns, each = length(df$years)),
                           model = 'linear')



# SSBtest <- seq(1, alpha/beta, length.out = 100)
# bhmodel <- (alpha*SSBtest)/(1+beta*SSBtest)
# 
# 
# 
# bhmodel <- alpha*SSBtest *exp(-beta * SSBtest)#*exp(-0.5*df$b[yr]*SDR^2+Ry)#
# plot(SSBtest,bhmodel, type = 'l')
# lines(rep(1/beta, 100), y = seq(0,1e12, length.out =100))



for(i in 1:nruns){
  set.seed(seeds[i])
  
  
  df <- load_data_seasons(nseason = 1,
                          nyear = 200,# Set up parameters 
                          Linf = 150, 
                          maxage = 10,
                          h = 0.4,
                          K = .4, 
                          t0 = 0, 
                          SDR = 0.2, # Recruitment deviations 
                          fishing.type = 'constant',
                          mortality = 'constant',
                          alpha = alpha,
                          beta = beta,
                          recruitment = 'Ricker',
                          negg = codest$parameters[['alpha.lin']],
                          eggbeta = codest$parameters[['beta.lin']],
                          F0 = 0) # Specify parameters
  df$egg.size <- df$egg.size/100
  
  
  
  tmprun <- run.agebased.true.catch(df, seed = seeds[i])
  
  
  df.save[df.save$run == i,]$SSB <- tmprun$SSB
  df.save[df.save$run == i,]$R <- tmprun$R.save
  df.save[df.save$run == i,]$Rtot <- tmprun$Rtot.save
  df.save[df.save$run == i,]$Catch <- tmprun$Catch
  
}


df.sum <- df.save %>% group_by(years, model) %>% 
  summarise(S = median(SSB),
            Rec = median(R),
            C = median(Catch),
            Rtot = median(Rtot),
            Smin = quantile(SSB, probs = 0.05),
            Smax = quantile(SSB, probs = 0.95),
            Recmin = quantile(R, probs = 0.05),
            Rmax = quantile(R, probs = 0.95),
            Cmin = quantile(Catch, probs = 0.05),
            Cmax = quantile(SSB, probs = 0.95)
            )



df.save.boff <- data.frame(years = rep(df$years, nruns),
                      SSB = NA,
                      R = NA,
                      Rtot = NA,
                      Catch = NA, 
                      run = rep(1:nruns, each = length(df$years)),
                      model = 'boff')

for(i in 1:nruns){
  set.seed(seeds[i])
  
  
  df <- load_data_seasons(nseason = 1,
                          nyear = 200,# Set up parameters 
                          Linf = 150, 
                          maxage = 10,
                          h = 0.4,
                          K = .4, 
                          t0 = 0, 
                          SDR = 0.2, # Recruitment deviations 
                          fishing.type = 'constant',
                          mortality = 'constant',
                          alpha = alpha,
                          beta = beta,
                          recruitment = 'Ricker',
                          negg = codest$parameters[['alpha.hyper']],
                          eggbeta = codest$parameters[['beta.hyper']],
                          F0 = 0) # Specify parameters
  df$egg.size <- df$egg.size/100
  
  
  tmprun <- run.agebased.true.catch(df, seed = seeds[i])
  
  df.save.boff[df.save.boff$run == i,]$SSB <- tmprun$SSB
  df.save.boff[df.save.boff$run == i,]$R <- tmprun$R.save
  df.save.boff[df.save.boff$run == i,]$Rtot <- tmprun$Rtot.save
  df.save.boff[df.save.boff$run == i,]$Catch <- tmprun$Catch
  
}


df.sum.boff <- df.save.boff %>% group_by(years, model) %>% 
  summarise(S = median(SSB),
            Rec = median(R),
            C = median(Catch),
            Rtot = median(Rtot),
            Smin = quantile(SSB, probs = 0.05),
            Smax = quantile(SSB, probs = 0.95),
            Recmin = quantile(R, probs = 0.05),
            Rmax = quantile(R, probs = 0.95),
            Cmin = quantile(Catch, probs = 0.05),
            Cmax = quantile(SSB, probs = 0.95)
  )




# Plot boff and not boff
df.plot <- rbind(df.save, df.save.boff)
df.sumplot <- rbind(df.sum, df.sum.boff)


ggplot(df.plot,aes(x = SSB, y = Rtot, color = model))+geom_point()+theme_classic()
  
  
ggplot(df.plot, aes(x= years, y = SSB, group = as.factor(run), color = model))+
  geom_line(alpha = 0.1, size = .5)+
  theme_bw()+theme(legend.position = 'none')+scale_y_log10()+
  geom_line(data = df.sumplot, aes(y = S, group = NA), size = 2)



ggplot(df.sumplot, aes(x = years, y = S, color = model))+geom_line()+theme_classic()+
  geom_hline(aes(yintercept = tmprun$SSB_0))+
  geom_ribbon(aes(ymin = Smin, ymax = Smax), fill = 'red', alpha = 0.1, linetype = 0)



plot(df.sum$Rtot/df.sum.boff$Rtot)
lines(df.sum.boff$Rtot)

plot(df.save$SSB, df.save.boff$SSB)

