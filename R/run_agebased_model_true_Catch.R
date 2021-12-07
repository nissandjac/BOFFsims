# Run an agebased model based on the parameters from the estimation model

run.agebased.true.catch <- function(df, seed = 123){
  # Constant parameters
  # Initialize the stock by age and sex
  
  # rSums <- function(x){
  #   if(df$>)
  #   
  set.seed(seed)  
  
  nseason <- df$nseason
  
  df$tEnd <- length(df$years)*nseason
  nyear <-df$tEnd/df$nseason # 
  year <- df$years
  tEnd <- nyear*nseason
  # True values 
  F0 <- df$F0
  
  # Set up the data spatially 
  nspace <- df$nspace
  recruitmat <- df$recruitmat
  if(df$move == FALSE){
    recruitmat <- 1
  }
  
  movemat <- df$movemat
  move.init <- df$move.init
  # M selectivity 
  Msel <- df$Msel # no difference between males and females
  M0 <- df$M0
  
  M <- t(replicate(Msel, n = nyear))
  M <- M[1:nyear,]*M0
  
  
  
  SDR <- exp(df$logSDR)
  b <- rep(1, nyear)
  # Survey selectivity 
   #surv.sel <- getSelec(df$age,df$parms$psel_surv, df$Smin_survey, df$Smax_survey) # Constant over time
  surv.sel <- df$sel # Same as fishing selectivity 
  # Catchability 
  q <- exp(df$logq) # Constant over time
  surv.sd <- exp(df$parms$logSDsurv) # Survey error
  
  # Maturity 
  Mat.sel <- df$Matsel # Fecundity
  h <- exp(df$parms$logh)
  
  # Age 
  nage <- df$nage
  age <- df$age
  
  
  R0 <- exp(df$parms$logRinit)
  

  if(df$recruitment == 'BH_steep'){
    R0 <- exp(df$parms$logRinit)
  }
  
  
  if(df$recruitment == 'BH'){
    R0 <- df$alpha/df$beta
  }
  
  if(df$recruitment == 'Ricker'){
    R0 <- 1/df$beta
  }
  
  Mage <- cumsum(M[1,])#c(0,cumsum(M[1,1:(nage-1)]))
  
  # Calculate N0 based on R0
  mage <- nage # Max age
  agetmp <- 0:(mage)
  nagetmp <- mage*1# SS multiplies this by 3 
  Mtmp <- rep(NA, nagetmp)
  
  N0 <- rep(NA,nagetmp)

  N0[1] <- R0
  N0[2:(nagetmp-1)] = R0*exp(-Mage[1:(nagetmp-2)])
  N0[nagetmp] =  R0*exp(-Mage[nagetmp-1])/(1-exp(-M[1,nage]))

  
  
  # Mtmp[1:nage] <- cumsum(M[1,])
  # Mtmp[nage+1:nagetmp] <- cumsum(M[1,])[nage]
  # 
  # N0tmp <- rep(NA,nagetmp)
  # 
  # N0tmp[1] <- R0
  # N0tmp[2:(nagetmp-1)] = R0*exp(-M[2:(nagetmp-1)]*age[2:(nagetmp-1)])
  # N0tmp[nagetmp] =  R0*exp(-M[nagetmp-1]*age[nagetmp-1])/(1-exp(-M[1,nage]))
  # 
  # N0 <- matrix(NA,nage)
  # N0[1:(nage-1)] <- N0tmp[1:(nage-1)]
  # N0[nage] <- sum(N0tmp[nage:nagetmp])
  # 
  SSB_0 <- NA
  Rprod0 <- NA
  
  for(i in 1:nspace){
    #SSB_0[i] <- sum(df$Matsel*N0*move.init[i])
    SSB_0[i] <- sum(N0*df$wage_ssb[,1]*df$Matsel)
    Rprod0[i] <- sum(N0*df$Matsel*df$egg.size)
  }
  names(SSB_0) <- paste(rep('space',each = df$nspace),1:nspace)
  
  R_0 <- R0#*move.init
  #R0_tot <- R0
  # Used the inital recruitment devs to get a start
  
  
  # Ninit <- rep(NA,nage)
  # Ninit_dev <- (df$parms$initN)
  # Ninit[1] <- R0
  # Ninit[2:(nage-1)] <-R0 * exp(-M[2:(nage-1)]*age[2:(nage-1)])*exp(-0.5*SDR^2*0+Ninit_dev[1:(nage-2)])
  # Ninit[nage] <- R0*exp(-(M[nage-1]*age[nage-1]))/(1-exp(-M[nage]))*exp(-0.5*SDR^2*0+Ninit_dev[nage-1])# Plus group (ignore recruitment dev's in first year )
  # 
  # Create containers to save the data
  # SSB_init <- NA
  # 
  # for(i in 1:nspace)
  # 
  
  #Ninit[1] <- sum((4*h*R_0*SSB_init/(SSB_0*(1-h)+ SSB_init*(5*h-1)))*exp(-0.5*1*SDR^2+df$parms$Rin[1]), na.rm = T)
  year_1 <- c(year,max(year)+1)

  SSB <- matrix(NA,nyear, nspace, 
                dimnames = list(year = df$years,
                               space = 1:nspace))
  SSB.all <- array(NA, dim = c(nyear, nspace,nseason),
                   dimnames = list(year = year, space = 1:nspace, season = 1:nseason))
  SSB.weight <- matrix(NA,nyear, nspace,
                       dimnames = list(year = year, space = 1:nspace))
  SSB.eq <- matrix(NA, nyear, nspace)
  
  
  Biomass.save <- matrix(NA,nyear, nspace, 
                         dimnames= list(year = year, space = 1:nspace))
  Catch <- matrix(NA,nyear, dimnames = list(year = year))
  Catch.age <- matrix(NA,nage,nyear, dimnames = list(age = age, year = year))
  CatchN <- matrix(NA,nyear, dimnames = list(year = year))
  CatchN.age <- matrix(NA,nage,nyear, dimnames = list(age =age, year = year))

  
  R.save <- matrix(NA,nyear, nspace, dimnames = list(year = year, space = 1:nspace))
  Rtot.save <- matrix(NA,nyear, nspace, dimnames = list(year = year, space = 1:nspace))
  
  Fsel.save <- array(NA,dim = c(nage,nyear,nspace), dimnames = list(age = age, year = year, space = 1:nspace))
  Fseason.save <- array(NA,dim = c(nage, nyear, nspace,nseason), dimnames = list(age = age, year = year, space = 1:nspace,
                                                                                 season = 1:nseason))
  Fout.save <- array(NA, dim = c(nyear,nseason,nspace), 
                     dimnames = list(year = year, season = 1:nseason, space = 1:nspace))
  
  N.save.age <- array(NA,dim = c(nage,nyear+1, nspace, nseason), 
                      dimnames = list(age = age, year = year_1, space = 1:nspace, season = 1:nseason))
  N.save.age.mid <- array(NA,dim = c(nage,nyear+1, nspace, nseason), 
                          dimnames = list(age = age, year = year_1, space = 1:nspace, season = 1:nseason))
  R.save <- matrix(NA, nyear, nspace)
  V.save <- array(NA,dim = c(nyear, nspace, nseason), dimnames = list(
    year = year, space = 1:nspace, season = 1:nseason))

  Catch.save.age <- array(NA,dim = c(nage,nyear, nspace, nseason), 
                          dimnames = list(age = age, year = year, space = 1:nspace, season = 1:nseason))
  CatchN.save.age <- array(NA,dim = c(nage,nyear, nspace, nseason), 
                           dimnames = list(age = age, year = year, space = 1:nspace, season =1:nseason))
  Catch.quota <- array(NA, dim = c(nyear, nspace, nseason), 
                       dimnames = list(year = year, space = 1:nspace, season =1:nseason))
  Catch.quota.N <- array(0, dim = c(nyear, nspace, nseason), dimnames = list(year = year, space = 1:nspace,
                                                                             season = 1:nseason))
    
  survey <- array(NA,dim = c(nyear), dimnames = list(year = year))
  survey.true <- array(NA, dim = c(nspace, nyear), dimnames = list(space = 1:nspace, year = year))
  surv.tot <- matrix(NA, nyear,nspace, dimnames = list(year = year, space = 1:nspace))

  age_comps_surv <- array(NA, dim = c(df$age_maxage,nyear), dimnames = list(age = 1:df$age_maxage,
                                                                             year = year)) # 
  age_comps_surv_space <- array(NA, dim = c(df$age_maxage,nyear,nspace), dimnames = list(
    age = 1:df$age_maxage, year = year))
 
  N.survey <- matrix(NA,df$age_maxage ,nyear, dimnames = list(age = 1:df$age_maxage,
                                                              year= year))

  age_comps_catch <- array(NA, dim = c(df$age_maxage,nyear), dimnames = list(age = 1:df$age_maxage,
                                                                             year = year))
  age_comps_catch_space <- array(NA, dim = c(df$age_maxage,nyear,nspace), dimnames = list(
    age = 1:df$age_maxage, year = year, space = 1:nspace))

  age_comps_OM <- array(NA, dim = c(nage,nyear, nspace,nseason), 
                        dimnames = list(age = age, year= year, space = 1:nspace, season = 1:nseason))
  
  Z.save <- array(NA, dim = c(df$nage, nyear,nspace,nseason), dimnames = list(age= age, year = year, space = 1:nspace,
                                                                              season = 1:nseason))

  Z.save[,1,1,1] <- M[1,]

  M.save <- array(NA, dim = c(df$nage, nyear,nspace,nseason), dimnames = list(age= age, year = year, space = 1:nspace,
                                                                              season = 1:nseason))
  M.save[,1,1,1] <- M[1,]
  
  Catch.age[,1] <- 0 # Assumed no fishing before data started 
  Catch[1]<- 0
  
  CatchN[1] <- 0
  CatchN.age[,1] <- 0
  
  survey[1] <- 1 # Surveys start later
  
  for (space in 1:nspace){
  survey.true[space,1] <- sum(N.save.age[,1,space,df$surveyseason]*surv.sel*q*df$wage_survey[,1])
  }

  idx.save <- seq(1,tEnd, by = nseason)
  
  # Distribute over space 
  Ninit <- rep(NA,nage)
  names(Ninit) <- age
  Ninit_dev <- (df$parms$initN)
  
  # Fishing before model is run
  Fpast <- df$Fpast*df$sel
  Fage_past <- cumsum(Fpast)
  
  if(all(is.na(Fage_past))){
    Fage_past <- rep(0, nage)
  }
  
  # if(recruitment != 'BH_steep'){
    Rprod0 <- R0
  # }
  
  
  
  Ninit[1] <- R0
  Ninit[2:(nage-1)] <-R0 * exp(-(Mage[1:(nage-2)]+Fage_past[1:(nage-2)]))*exp(-0.5*SDR^2*0+Ninit_dev[1:(nage-2)])
  Ninit[nage] <- R0*exp(-((Mage[nage-1]+Fage_past[nage-1])))/(1-exp(-(M[1,nage]+Fpast[nage])))*exp(-0.5*SDR^2*0+Ninit_dev[nage-1])# Plus group (ignore recruitment dev's in first year )
  
  #p.save <-matrix(NA,tEnd)
  # SSB_0[i] <- sum(Ninit*df$wage_ssb[,1]*df$Matsel)
  # Rprod0[i] <- sum(N0*df$Matsel*df$egg.size)
  # 
  
  
  
  for (space in 1:nspace){
      # if (season == 1){
      N.save.age[,1,space,1] <- Ninit#*move.init[space] # Just to initialize 
      N.save.age.mid[,1,space,1] <- N.save.age[,1,space,1]*exp(-0.5*(M[1,]/nseason))
      # }else{
      #   N.save.age[,1,space,season] <- N.save.age[,1,space,season-1]*exp(-M/nseason)
      #   N.save.age.mid[,1,space,season] <- N.save.age[,1,space,season]*exp(-0.5*(M/nseason))
      # }
      # }
  }
  
  
  
  
  
  #Fspace <- c(0.2612,0.7388) # Contribution of Total catch (add to one)    #Z <- (Fyear+Myear)
  Fnseason <- df$Fnseason
  
  if(sum(Fnseason) != 1){
    warning('Catch input does not match')
  }
  
  pope.mul <- nseason/1*0.5
  pope.mul <- 0.50
  
  # Recruitment paramerters 
  beta <- df$beta
  recruitment <- df$recruitment
  
  
  
 
  for (yr in 1:nyear){ # Loop over years add one year for initial distribution
    
    w_catch <- df$wage_catch[,yr]
    w_surv <- df$wage_survey[,yr]
    w_mid <- df$wage_mid[,yr]
    w_ssb <- df$wage_ssb[,yr]
  
    
    
    Ry <- df$parms$Rin[yr]
 

    # Fyear <- F0[yr]*Fsel
    Myear <- M0[yr]
    
    if(df$move == FALSE){
      Fspace <- 1 # All catches 
    }
    
    Mseason <- Myear/nseason # M is distributed throughout the year
    
    
    # fix Ssb and recruitment in all areas 
    for(space in 1:nspace){
      SSB.weight[yr,space] <- sum(N.save.age[,yr,space,1]*as.numeric(w_ssb)*df$Matsel, na.rm = TRUE)
      SSB[yr,space] <- sum(N.save.age[,yr,space,1]*as.numeric(w_ssb)*df$Matsel, na.rm = TRUE)
      Biomass.save[yr,space] <- sum(N.save.age[,yr,space,1]*as.numeric(w_mid), na.rm = TRUE)
      SSB.all[1,space,1]<- sum(N.save.age[,1,space,1]*Mat.sel, na.rm = TRUE)
      
      # Calculate the unfished spawning biomass given the annual life history params
      Mage <- cumsum(Myear)#c(0,cumsum(M[1,1:(nage-1)]))
      
      # Calculate N0 based on R0
      # N0.eq <- rep(NA,nagetmp)
      # 
      # N0.eq[1] <- R0
      # N0.eq[2:(nagetmp-1)] = R0*exp(-Mage[1:(nagetmp-2)])
      # N0.eq[nagetmp] =  R0*exp(-Mage[nagetmp-1])/(1-exp(-Myear[nage]))
      # SSB.eq[yr, space] <- sum(N0.eq*as.numeric(w_ssb)*df$Matsel)
      # 
    # Recruitment only in season 1  
      if(recruitment == 'BH'){
      
      N.save.age[1,yr,space,1] <- 0 # No recruits yet (avoid NA calculation)
      # Calculate eggs per individual   
      Rtot <- sum(N.save.age[,yr,space,1]*Mat.sel*df$egg.size)
        
      # R <- (4*h*R_0[space]*SSB[yr,space]/
      #       (SSB_0[space]*(1-h)+ SSB[yr,space]*(5*h-1)))*exp(-0.5*df$b[yr]*SDR^2+Ry)#*recruitmat[space]
    
      # Formulate BH in terms of egg production instead of SSB 
      R <- (df$alpha*Rtot)/(1+df$beta*Rtot)*exp(-0.5*df$b[yr]*SDR^2+Ry)
      
   
      
      }
      
      if(recruitment == 'BH_steep'){
        
        N.save.age[1,yr,space,1] <- 0 # No recruits yet (avoid NA calculation)
        # Calculate eggs per individual   
        Rtot <- sum(N.save.age[,yr,space,1]*Mat.sel*df$egg.size)

        # 
        R <- (4*h*R_0*Rtot/
               (R_0*(1-h)+ Rtot*(5*h-1)))*exp(-0.5*df$b[yr]*SDR^2+Ry)#*recruitmat[space]

        # 
        R <- (4*h*R_0[space]*SSB[yr,space]/
                (SSB_0[space]*(1-h)+ SSB[yr,space]*(5*h-1)))*exp(-0.5*df$b[yr]*SDR^2+Ry)#*recruitmat[space]

        print('hm fix issues with this SR function')
      }
      
      
      if(recruitment == 'Ricker'){

        N.save.age[1,yr,space,1] <- 0 # No recruits yet (avoid NA calculation)
        # Calculate eggs per individual   
        Rtot <- sum(N.save.age[,yr,space,1]*Mat.sel*df$egg.size)
        #Rtot <- SSB[yr,space]
        
        R <- alpha*Rtot *exp(-beta * Rtot)*exp(-0.5*df$b[yr]*SDR^2+Ry)#
     
        
        
        
      }
      
      if(recruitment == 'linear'){
        
        Rtot <- sum(N.save.age[,yr,space,1]*Mat.sel*df$egg.size, na.rm = TRUE)
        #Rtot <- SSB[yr,space]
        
        R <- Rtot*exp(-0.5*df$b[yr]*SDR^2+Ry)#
        
      }
      
      
      if(recruitment == 'BH_R'){
        
        
        Rtot <- sum(N.save.age[,yr,space,1]*Mat.sel*df$egg.size, na.rm = TRUE)
        
        
        # Regulate R0 based on the lambda function 
        x <- sum(N.save.age[df$age >= df$tau,yr,,]*df$Matsel[df$age >= df$tau]*df$egg.size[df$age >= df$tau])/
          sum(N.save.age[,yr,,]*df$Matsel*df$egg.size, na.rm  =TRUE)
        
        
        if(is.na(df$lambda)){
          R0_boff <- R0
        }else{
         R0_boff <- R0*1/(exp(-df$lambda.slope*(x-lambda)))
        }
        
        
        
        
        # if(yr == 1){
        #   print(df$lambda.slope)
        # }
        
        # print(round(R0_boff/R0, 3))

        R <- R0_boff*Rtot/(Rtot+R0_boff)*exp(-0.5*df$b[yr]*SDR^2+Ry)#
        
        
      }
      
      
      
      if(R < 1e-10){ # Remove numerical instability
        R <- 0
      }
      
      
      N.save.age[1,yr,space,1] <- R
      R.save[yr,space] <- R
      Rtot.save[yr,space] <- Rtot
      
      
      
    }
    
  
    
    
    for (season in 1:nseason){
      for (space in 1:nspace){
        

       
        #p.save[yr] <- sum(pseltmp)
        # 
        Fsel <- df$sel

        Fsel.save[,yr,space] <- Fsel
        
        if(nspace > 1){
          if(df$years[yr]<= 2018){
            Catch_space <- df$Catch.country[yr,space]
          }else{
            Catch_space <- df$Catch[yr]*Fspace[space]  
          }
        }else{
          Catch_space <- df$Catch[yr]
        }
        
        
        B.tmp <-  sum(N.save.age[,yr,space,season]*exp(-Mseason*pope.mul)*w_catch*Fsel) # Get biomass from previous year
        N.tmp <- N.save.age[,yr,space,season]#
        V.save[yr,space,season] <- B.tmp
        
        if(df$catch.true == TRUE){
          E.temp <- Catch_space*Fnseason[season,space]#*Fspace[space] # Catch distribution in the year
          Catch.quota[yr,space,season] <- E.temp
        
        
        if(E.temp/B.tmp >= .95){
          
          warning(paste('Catch exceeds available biomass in year:',year[yr],' and season', season, 'area', space)) # Stop if in the past 
                    #print(paste('Catch exceeds available biomass in year:',year,' and season', season, 'area', space))
          E.temp <- 0.95*B.tmp
          Catch.quota.N[yr,space,season] <- 1
          #if(df$years[yr] > 2026){
          #stop('danger')
          #  }
          
        }
        
        Fout <- getF(E.temp,B.tmp,Mseason = Mseason, Fsel = Fsel, N.tmp = N.tmp, w_catch = w_catch, 
                     method = 'Hybrid')
        
        if(E.temp>0){
          Fseason <- Fout*Fsel
          Fnew <- Fout  
          Z <- Fnew*Fsel+Mseason
          Fseason <- Fnew*Fsel
        }else{
          Fseason <- 0
        }
        
        }else{
          Fout = df$F0[yr]/nseason
          Fnew <- Fout  
          Z <- Fnew*Fsel+Mseason
          Fseason <- Fnew*Fsel
        }
        
        #Fout <- df$parms$F0[yr]
        
      
        Fout.save[yr,season,space] <- Fout # terminal fishing mortality 
        
        Fseason.save[,yr,space,season] <- Fseason
      
        Z <- Mseason+Fseason
        Z.save[,yr,space,season]<- Z
        M.save[,yr,space,season] <- Mseason
        
        if(((space-1) == 0)){
          spaceidx <- 2
        }
        if(space == nspace){
          spaceidx <- nspace-1
        }
        if(space > 1 & space < nspace){
          spaceidx <- c(space-1,space+1)
        }
        
        if(df$move == FALSE){
          spaceidx <- 1
        }
        
        if(season <nseason){
        
          N.save.age[,yr,space,season+1] <- N.save.age[,yr,space,season]*exp(-Z)-
          N.save.age[, yr,space,season]*exp(-Z)*(movemat[space,,season,yr])+ # Remove the ones that leave
          N.save.age[, yr,spaceidx,season]*exp(-Z)*(movemat[spaceidx,,season,yr])# add the ones come to the surrounding areas
        
        age_comps_OM[,yr,space,season] <- N.save.age[, yr,space,season]/sum(N.save.age[, yr,space,season])
        
        SSB.all[yr,space,season]<- sum(N.save.age[,yr,space,season]*Mat.sel, na.rm = T)
        
        Catch.save.age[, yr,space, season] <- (Fseason/(Z))*(1-exp(-(Z)))*N.save.age[,yr,space,season]*w_catch
        CatchN.save.age[, yr,space, season] <- (Fseason/(Z))*(1-exp(-(Z)))*N.save.age[,yr,space,season]
        
        
        }else{
          N.save.age[2:(nage-1),yr+1,space,1] <- N.save.age[1:(nage-2),yr,space,season]*exp(-Z[1:(nage-2)])-
            N.save.age[1:(nage-2), yr,space,season]*exp(-Z[1:(nage-2)])*(movemat[space,1:(nage-2),season,yr])+ # Remove the ones that leave
            N.save.age[1:(nage-2), yr,spaceidx,season]*exp(-Z[1:(nage-2)])*(movemat[spaceidx,1:(nage-2),season,yr])# add the ones come to the surrounding areas
          
          # Plus group 
          Nsurvive.plus <- (N.save.age[nage-1, yr,space, nseason]*exp(-Z[nage-1])+
                              N.save.age[nage, yr,space, nseason]*exp(-Z[nage]))
          
          Nout.plus <- Nsurvive.plus*(movemat[space,nage, season,yr]) # Leaving

          
          Nin.plus <- (N.save.age[nage-1, yr,spaceidx,nseason]*exp(-Z[nage-1])+
                         N.save.age[nage, yr,spaceidx,nseason]*exp(-Z[nage]))*
                          (movemat[spaceidx,nage, season,yr]) # Incoming

          N.save.age[nage,yr+1,space,1] <- Nsurvive.plus- Nout.plus + Nin.plus
          
          
          age_comps_OM[,yr,space,season] <- N.save.age[, yr,space,season]/sum(N.save.age[, yr,space,season])
          
          SSB.all[yr,space,season]<- sum(N.save.age[,yr,space,season]*Mat.sel, na.rm = T)
          Catch.save.age[, yr,space, season] <- (Fseason/(Z))*(1-exp(-(Z)))*N.save.age[,yr,space,season]*w_catch
          CatchN.save.age[, yr,space, season] <- (Fseason/(Z))*(1-exp(-(Z)))*N.save.age[,yr,space,season]
         
          
        }
   
        
        if(is.na(SSB[yr,space])){
          stop('SSB is NA')
        }
      } 
      
      if(df$catch.true == 1){
      if(Catch.quota[yr,space,season]>0){
        if((sum(Catch.save.age[, yr,space, season])/Catch.quota[yr,space,season]) > 1.05){
         # stop('F estimation overshoots more than 10%')
        }
      }
      }
      } # End of season loop
    
  
    
    #Catch.age[,idx]  <- (Fyear/(Fyear+Myear))*(1-exp(-(Fyear+Myear)))*rowSums(N.save.age[,idx,,1])*w_catch # Calculate the catch in kg 
    
    if(nseason>1){
    Catch.age[,yr] <- apply(Catch.save.age[,yr,,],MARGIN = 1,FUN = sum)
    Catch[yr] <- sum(Catch.save.age[,yr,,])  
  
    CatchN.age[,yr] <- apply(CatchN.save.age[,yr,,],MARGIN = 1,FUN = sum)
    CatchN[yr] <- sum(CatchN.save.age[,yr,,])  
    }else{
      
      if(nspace == 1){
      Catch.age[,yr] <- Catch.save.age[,yr,,]
      Catch[yr] <- sum(Catch.save.age[,yr,,])
      
      CatchN.age[,yr] <- CatchN.save.age[,yr,,]
      CatchN[yr] <- sum(CatchN.save.age[,yr,,])
      }else{
        Catch.age[,yr] <- rowSums(Catch.save.age[,yr,,])
        Catch[yr] <- sum(Catch.save.age[,yr,,])
        
        CatchN.age[,yr] <- rowSums(CatchN.save.age[,yr,,])
        CatchN[yr] <- sum(CatchN.save.age[,yr,,])
      }
    }  
    
    Msurveymul <- df$smul # maybe later fix this 
        
    for (space in 1:nspace){
    survey.true[space,yr] <- sum(N.save.age[,yr,space,df$surveyseason]*
                                    exp(-Msurveymul*Z.save[,yr,space,df$surveyseason])*surv.sel*q*w_surv)
      
    }
    
    
    #  }
      # Save the survey 
      # Survey is conducted in the start of the year
    # }else{
    #   Msurveymul <- 0.5
    # }
    
    if(df$move == FALSE){
      Nsurv <- N.save.age[,yr,,df$surveyseason]*
        exp(-Msurveymul*Z.save[,yr,space,df$surveyseason])
    }else{
      Nsurv <- rowSums(N.save.age[,yr,,df$surveyseason]*
                         exp(-Msurveymul*Z.save[,yr,space,df$surveyseason]))
    }
    
    if (df$flag_survey[yr] == 1){
        
        err <- rnorm(n = 1,mean = 0, sd = surv.sd)
        surv <- exp(log(sum(Nsurv*surv.sel*q*w_surv))+err) # If the xtra factor is not included the mean is > 1
        survey[yr] <- surv
      }else{
        survey[yr] <- 1
      }
    
      Ntot.year <- Nsurv
      
      surv.tmp <- sum(Ntot.year*surv.sel*q)
      
     if(df$flag_survey[yr] == 1){
        age_comps_surv[1,yr] <- 0 # No year 1 recorded

        age_comps_surv[1:(df$age_maxage-1),yr] <-  (Ntot.year[2:(df$age_maxage)]*surv.sel[2:(df$age_maxage)]*q)/surv.tmp
        age_comps_surv[df$age_maxage,yr] <- sum(Ntot.year[(df$age_maxage+1):nage]*surv.sel[(df$age_maxage+1):nage]*q)/surv.tmp
      }else{
        age_comps_surv[,yr] <- NA
      }
      
      for(space in 1:nspace){
      Ntot.year <- N.save.age[,yr,space,df$surveyseason]
      surv.tot[yr,space]  <- sum(Ntot.year*surv.sel*q*exp(-Msurveymul*Z.save[,yr,space,df$surveyseason]))
        
      age_comps_surv_space[1,yr,space] <- 0 # No year 1 recorded
        
      age_comps_surv_space[1:(df$age_maxage-1),yr,space] <-  
        (Ntot.year[2:(df$age_maxage)]*surv.sel[2:(df$age_maxage)]*q)/surv.tot[yr,space]
      age_comps_surv_space[df$age_maxage,yr,space] <- 
        sum(Ntot.year[(df$age_maxage+1):nage]*surv.sel[(df$age_maxage+1):nage]*q)/surv.tot[yr,space]
      
      if(nseason>1){
      Catch.tmp <- rowSums(CatchN.save.age[, yr,space,])
      }else{
      Catch.tmp <- CatchN.save.age[, yr,space,]
      }
      
      Catch.tot <- sum(CatchN.save.age[,yr,space,])
      
      age_comps_catch_space[1:(df$age_maxage-1),yr,space] <- Catch.tmp[2:(df$age_maxage)]/Catch.tot
      age_comps_catch_space[df$age_maxage,yr,space] <- sum(Catch.tmp[(df$age_maxage+1):nage])/Catch.tot
      
      
      
      }
      
      #
      if(df$flag_catch[yr] == 1){
        age_comps_catch[1:(df$age_maxage-1),yr] <-  CatchN.age[2:(df$age_maxage),yr]/CatchN[yr]
        age_comps_catch[df$age_maxage,yr] <- sum(CatchN.age[(df$age_maxage+1):nage,yr])/CatchN[yr]
        
      }else{
        age_comps_catch[,yr] <- NA
      }
      

    }# End of year loop
  #}
  
    if(df$move == FALSE){
      Nsave <- N.save.age[,,,nspace]
      SSB.save <- SSB
    }else{  
    Nsave <- apply(N.save.age[,,,1],2,rowSums)
    SSB.save <- rowSums(SSB)
    }
  
    # Add names to output 
    year_1 <- c(df$years,max(df$years+1))

    # 

    

    df.out   <- list(N.save = Nsave, 
                     SSB = SSB, 
                     N.save.age = N.save.age,
                     R.save = R.save,
                     Rtot.save = Rtot.save,
                     V.save = V.save,
                     Biomass = Biomass.save,
                     SSB.all = SSB.all,
                     SSB.eq = SSB.eq,
                     Catch.save.age = Catch.save.age,
                     CatchN.save.age = CatchN.save.age,
                     Catch = Catch, 
                     Catch.age = Catch.age, 
                     Catch.quota = Catch.quota,
                     Catch.quota.N = Catch.quota.N,
                     Fout = Fout.save,
                     age_comps_OM = age_comps_OM,
                     age_catch = age_comps_catch,
                     SSB_0 = SSB_0, 
                     N0 = N0,
                     SSB.weight = SSB.weight,
                     survey.true = survey.true,
                     Z = Z.save,
                     M = M.save,
                     survey = as.numeric(survey),
                     age_comps_surv = age_comps_surv,
                     age_comps_country = age_comps_surv_space,
                     age_comps_catch_space = age_comps_catch_space,
                     Fseason = Fseason.save,
                     Fsel = Fsel.save, 
                     Ninit = Ninit,
                     R0 = Rprod0)

  return(df.out)
    
}

