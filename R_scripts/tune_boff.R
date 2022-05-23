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

F0 <- 0.1
nruns <- 30
years <- 100


rho <- c(.1,.9)#, 0.3, 0.5, .9)
set.seed(123)

lambda <- 3


df <- load_data_seasons(nseason = 1,
                        nyear = years,# Set up parameters 
                        Linf = Linf, 
                        maxage = maxage,
                        K = K, 
                        h = 1,
                        t0 = 0, 
                        M = M0,
                        tau_sel = tau_sel,
                        tau = tau,
                        SDR = 0, # Recruitment deviations 
                        fishing.type = fishing.type,
                        mortality = mortality,
                        alpha = alpha,
                        beta = beta,
                        recruitment = recruitment,
                        recruitment.type = recruitment.type,
                        F0 = seq(0,.9, length.out = years),
                        SDF = .7,
                        rhoR = rho[1],
                        R0 = R0,
                        lambda = lambda,
                        lambda.cut = .9) # Specify parameters






# Do boff as a function of autocorrelation 
boff <- run.agebased.true.catch(df, seed = 123)

plot(tmp$x.save, tmp$R.save)
df$lambda <- NA

noboff <- run.agebased.true.catch(df, seed = 123)
plot(noboff$x.save, noboff$R.save, ylim = c(500, 1500))
lines(boff$x.save, boff$R.save)


df.plot <- data.frame(yr = rep(df$years,2),  frac = c(noboff$x.save, boff$x.save),
                      R = c(noboff$R.save, boff$R.save),
                      F0 = c(df$F0, df$F0),
                      Catch = c(noboff$Catch, boff$Catch),
                      Rtot = c(noboff$Rtot.save, boff$Rtot.save),
                      SSB = c(noboff$SSB, boff$SSB),
                      model = rep(c('noBoff','BOFF'), each = df$nyear)
)


p1 <- ggplot(df.plot, aes(x = frac, y = R, color = model))+geom_line(size = 1.3)+theme_classic()+
  theme(legend.position = c(0.2,.8))+scale_x_continuous('fraction of spawners which age are over 90% mature')

p1


ggsave(filename = 'fractionSpawner.png', plot = p1, width = 16, height = 10, units = 'cm')




# Compare Fmsy and MSY 
p2 <- ggplot(df.plot, aes(x = F0, y = Catch, color = model))+geom_line()+theme_classic()
p3 <- ggplot(df.plot, aes(x = SSB, y = Catch, color = model))+geom_line()+theme_classic()

p2/p3


ggsave(filename = 'fishingstuff.png', plot =p2/p3, width = 16, height = 10, units = 'cm')



### Run a couple of scenarios 
models <- 'linear'
recLambda <- c('BOFF', 'noBOFF')
F0 <- 0.1
fishing.type <- 'AR'
Fpast <- F0
lambda.in <- .5
rho <- c(0.1, .9)
t0 <- 0
lambda.cut <- 0.8
SDR <- c(.2,.6,1.2)
SDF = .2
M <- .3
nruns

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
                        SDR = SDR,
                        F0 =F0 ,
                        SDF = SDF,
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




df.plot <- ls.plot[[2]]


ggplot(df.plot[df.plot$years > 50,] ,aes(x = SSB, y = R, color = model))+geom_point()+facet_wrap(~rho)+
  theme_classic()

ggplot(df.plot[df.plot$years > 20,] ,aes(x = xfrac, y = R, color = model))+geom_point(alpha = 0.2)+facet_wrap(~rho)+
  theme_classic()+geom_smooth(method = 'lm', size = 2)+scale_y_log10()

df.p2 <- df.plot[df.plot$years > 20,] %>%  group_by(model, run, rho) %>% summarise(SSBcor = cor(SSB, R),
                                                                              xCor = cor(xfrac, R))

p1 <- ggplot(df.p2, aes(x = model, y = SSBcor, color = model))+geom_violin()+geom_boxplot(width = .2)+facet_wrap(~rho)+theme_bw()
p2 <- ggplot(df.p2, aes(x = model, y = xCor, color = model))+geom_violin()+geom_boxplot(width = .2)+facet_wrap(~rho)+theme_bw()


p1 / p2

