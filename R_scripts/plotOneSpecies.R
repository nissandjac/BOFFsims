## PLot one species ## 

load('all_runs.Rdata')


# Make sure we dont do anything weird # 

df.save <- ls.plot[[2]]

p1 <- ggplot(df.save, aes(x = xfrac, y = R, color = model))+facet_wrap(~rho)+
  geom_point(alpha = 0.05)+geom_smooth()
p1



# Caclulate the fraction of old invividuals 
# 
# df.N <- ls.plot[[3]]
# df.N$old <- 'old'
# df.N$old[df.N$age < 5] <- 'young'
# df.N[is.na(df.N$lambda),]$lambda <- 0
# 
# df.Ntot <-  df.N %>% group_by(years, rho, lambda, run, model, old) %>%
#   summarise(SSB = sum(N*weight*mat),
#             F0 = mean(F0),
#             Catch = sum(Catch)
#   ) %>% arrange(lambda, rho,run , years, model)
# 
# 
# young <- df.Ntot[df.Ntot$old == 'young',]
# old <- df.Ntot[df.Ntot$old == 'old',]
# 
# tot <- data.frame(years = old$years,
#                   rho = old$rho,
#                   run = old$run,
#                   model = old$model,
#                   lambda = old$lambda,
#                   oldFrac = old$SSB/(old$SSB+young$SSB),
#                   Catch = old$Catch+young$Catch)
# 
# # Now add recruitment and devs 
# 
# R <- ls.plot[[2]]
# 
# # Check that the data frames are ordered correctly 
# plot(tot$Catch[1:2000])
# lines(R$Catch[1:2000])
# 
# range(tot$Catch[1:300]/R$Catch[1:300])
# 
# head(tot[101:250,])
# head(R[101:250,])
# 
# 
# # now combine 
# tot$xfrac <- R$xfrac
# tot$R <- R$R
# tot$Rdev <- R$Rdev
# tot$Rscam <- R$Rscam
# tot$Rscam_dev <- R$Rscam_dev
# 
# # Now do the analysis for violin plots again 
# 
tot.cor <- df.save %>% group_by(rho, run, model) %>% summarise(R_old = cor(R, xfrac),
                                                           Rdev_old = cor(Rdev, xfrac),
                                                           scam_old = cor(Rscam, xfrac),
                                                           scamdev_old = cor(Rscam_dev, xfrac)) %>%
  pivot_longer(c(R_old, Rdev_old, scam_old, scamdev_old))
# 

dodge <- position_dodge(.5)


p1 <- ggplot(tot.cor, aes(x = name, y = value, color = model, fill = model))+geom_violin(position = dodge, alpha = .3)+
  geom_boxplot(width = 0.2,position = dodge, col ='black', show.legend = F)+theme_classic()+geom_hline(aes(yintercept = 0), linetype = 2)+
  facet_wrap(~rho)+scale_x_discrete('', labels = c('cor\n(R, x)',
                                                   'cor\n(Rdev, x)',
                                                   'cor\n(R (scam), x)',
                                                   'cor\n(Rdev (scam), x)'))+
  theme(legend.position = 'top',
        axis.text.x = element_text(angle = 35, vjust = .5, hjust=0),
        legend.title = element_blank())
p1

ggsave(file = 'correlation_plot_recruitment.png', p1, width = 16, height = 10, units = 'cm')













  

