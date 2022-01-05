### Plot the four different recruitment functions ### 
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


# Punch in life history parameters (these are standard in run scenario)

df <- load_data_seasons()

# BH first 


h <- .8
R0 <- exp(df$parms$logRinit)
alpha <- 10
beta <- alpha/R0
SSB0 <- calcSSB0(df)

Rtot <- seq(1,R0*2, length.out = 150) 
R_bh<- (alpha*Rtot)/(1+beta*Rtot)
  
plot(Rtot, R_bh, type = 'l', ylim = range(c(Rtot,R_bh)))
lines(Rtot,rep(R0, 150))

R_steep <- (4*h*R0*Rtot)/
          (R0*(1-h)+ Rtot*(5*h-1))
  
plot(Rtot, R_steep, type = 'l')

beta <- 3e-4
alpha <- 7
R_ricker <- alpha*Rtot *exp(-beta * Rtot)
  
plot(Rtot, R_ricker, type='l')

R_eggs <- R0*Rtot/(Rtot+R0)

  
plot(Rtot,R_eggs, type = 'l')  

df.plot <- data.frame(Eggs = rep(Rtot, 4),
                      Recruits = c(R_bh, R_steep, R_ricker, R_eggs),
                      R0 = R0,
                      model = rep(c('Beverton Holt','Steepness' , 'Ricker', 'Egg BH'), each = 150))   
df.plot$R0[df.plot$model == 'Ricker'] <- NA


p1 <- ggplot(df.plot, aes(x= Eggs, y = Recruits, color = model))+geom_line()+facet_wrap(~model, nrow = 4)+theme_classic()+
  geom_line(aes(y = R0), col = 'black', linetype = 2)+theme(legend.position = 'none')

p1

## Plot different autocorrelation in recruitment 

rhoR <- c(0.1, 0.4, 0.7, 0.9)
SDRin <- c(0.1, 0.5, 0.8, 1)

  rhodf <- data.frame(rho = rep(rhoR, each = df$nyear*4),
                      SDR = rep(SDRin, each = df$nyear),
                 recdev = NA, 
                 time = 1:df$nyear)

for(i in 1:length(rhoR)){
  for(j in 1:length(SDRin)){
    tmp <- load_data_seasons(rhoR = rhoR[i], 
                           SDR = SDRin[j], 
                           recruitment.type = 'AR')
  
  
  rhodf[rhodf$rho == rhoR[i] & rhodf$SDR == SDRin[j],]$recdev <- (tmp$parms$Rin)
  }
}

p2 <- ggplot(rhodf, aes(x = time, y = (recdev), color = as.factor(rho)))+geom_line(size = 1.2)+facet_grid(SDR~rho)+theme_classic()+
  scale_y_continuous('Recruitment deviations')+coord_cartesian()+ggtitle('Autocorrelation coefficient')+
  theme(legend.position = 'none')



png('figures/rec_types.png', res = 400, width = 10, height = 20, units = 'cm')
p1
dev.off()


png('figures/AR_recruitment.png', res = 400, width = 16, height = 16, units = 'cm')
p2
dev.off()


## Illustrate Lambda 
x <- seq(0.01, .99, length.out = 50)
lambda <- 0.4
lambda.slope <- 0.7


R0 <- 100*1/(exp(-lambda.slope*(x-lambda)))

plot(x, R0)

df.plot <- data.frame(x = x, R0 = R0)


png('figures/lambdacalc.png', width = 8, height = 8, res = 400, units = 'cm')
ggplot(df.plot, aes(x = x, y = R0))+geom_line(size = 1.2)+theme_classic()+geom_hline(aes(yintercept = 100), linetype = 2)+
  geom_vline(aes(xintercept = lambda), linetype = 2)
dev.off()




