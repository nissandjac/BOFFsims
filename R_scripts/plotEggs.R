# Plot the diferent eggs 

df <- read.csv('data/fecundityEggSizeFemaleSize.csv')



ggplot(df, aes(x = FemaleSize_mm, y = Fecundity_nOfEggs_per_female, color = Species))+
  geom_point()+scale_y_continuous()+scale_x_log10()+facet_wrap(~Species, scales = 'free')+
  geom_smooth(method = 'lm', aes(group = Species))+theme_bw()+theme(legend.position = 'none')


# Just use cod as an example 

cod <- df[df$Species=="Gadus morhua" & df$Location=="Coastal Iceland",]


# Plot cod # 

ggplot(cod, aes(y = Fecundity_nOfEggs_per_female, x= FemaleSize_mm))+geom_point()


# Fit an exponential and a linear model # 

df.tmb <- list(
  w = (.01*(.1*cod$FemaleSize_mm)^3)/1000,
  fecundity = cod$Fecundity_nOfEggs_per_female,
  nobs = length(cod$Fecundity_nOfEggs_per_female)
)


parms <- list(lalpha = 60)

compile("src/fits.cpp")
dyn.load(dynlib("src/fits"))
obj <-MakeADFun(df.tmb,parms,DLL="fits")


plot(df.tmb$w,df.tmb$fecundity/1e6, xlim = c(0,30))
lines(df.tmb$w, df.tmb$w*1)


lower <- obj$par-Inf
upper <- obj$par+Inf

system.time(opt<-nlminb(obj$par,obj$fn,obj$gr,lower=lower,upper=upper, 
                        control = list(iter.max = 1e6,
                                       eval.max = 1e6))) #

