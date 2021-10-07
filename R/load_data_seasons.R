## Load the hake data
# year and age input 
load_data_seasons <- function(nseason = 1, 
                              nspace = 1, 
                              nyear = 50,
                              Linf = 30, 
                              K = 1, 
                              t0 = -0.1,
                              gamma = 1.1,
                              tau = 3,
                              gamma_sel = 1.1,
                              tau_sel = 3,
                              maxage = 10,
                              M = 0.2,
                              b = 1,
                              R0 = 1e4,
                              h = 0.8,
                              q = 1,
                              smul = 0.5,
                              SDsurv = 0.1,
                              SDcatch = 0.01,
                              movemaxinit = 0.35, 
                              movefiftyinit = 6,
                              nsurvey = 2, 
                              SDR = .6, 
                              Fpast = 0,
                              bfuture = 0.5,
                              moveout = 0.85, 
                              movesouth = 0.05,
                              moveinit = NA, 
                              moveslope = 0.9,
                              selectivity_change = 0,
                              yr_future  = 0,
                              sel_hist = 1,
                              catch.true = FALSE,
                              catch = NA,
                              F0 = 0.2,
                              fishing.type = 'constant',
                              SDF = 0.05,
                              rho_F = 0.7,
                              SDM = 0.1,
                              rho = 0.8,
                              M_limit = NA,
                              mortality = 'constant',
                              recruitment = 'BH',
                              alpha = NA,
                              beta = 1, 
                              negg = rep(1, nage),
                              eggbeta = 1,
                              seed = NA
                              ){
  
  #' @nseason = number of seasons 
  #' @nspace = Spatial areas
  #' @nyear = Final year that contains empirical data
  #' @Linf  Asymptotic length
  #' @K VonB growth rate
  #' @t0 VonB age at 0
  #' @gamma slope of maturity function
  #' @tau age at 50% maturity 
  #' @gamma_sel slope of selectivity function 
  #' @tau_sel age at 50% selectivity 
  #' @maxage maximum age of stock
  #' @M vector or integer of natural mortality
  #' @b bias adjustment (scalar or vector)
  #' @R0 unfished recruitment 
  #' @h steepness
  #' @q survey catchability
  #' @smul Annual survey timing (i.e., 0.5 middle of the year)
  #' @SDsurv survey uncertainty 
  #' @SDcatch catch uncertainty 
  #' @movemaxinit = max movement rate 
  #' @movefiftyinit = age at 50% movement rate
  #' @nsurvey survey frequency, i.e., 2 = every second year
  #' @SDR Standard deviation of recruitment
  #' @Fpast Fishing mortality before beginning of run 
  #' @bfuture recruitment bias adjustment in the future - scalar
  #' @moveout fraction of individuals that travel south in the last season
  #' @movesouth fraction of individuals that move south during the year
  #' @moveinit Initial distribution of fish
  #' @moveslope Slope of the movement function
  #' @selectivity_change flag for selectivity changing in the future 0
  #' @yr_future Create dummy data for future years
  #' @sel_hist use historical hake selecetivity 1 or 0 value
  #' @catch.true Use fishing mortality or catch?
  #' @F0 Fishing mortality input vector or integer, if catch is not provided
  #' @fishing.type how is fishing mortality modeled (constant or AR)  
  #' @SDF Variation in AR fishing mortality 
  #' @rho_F autocorrelation in fishing mortality 
  #' @SDM natural mortality random walk sd 
  #' @rho natural mortality autocorrelation parameter
  #' @M.type which type of natural mortality 1) constant, 2) random walk, 3) increasing, 4) decreasing

  if(is.na(seed)){
    set.seed(round(runif(1, min = 1, max = 1e6)))
  }else{
    set.seed(seed)
  }
  
  if(is.na(M_limit)){
    M_limit <- 2*M
  }
  
  years <- 1:nyear
  nyear <- length(years)
  tEnd <- length(years)*nseason
  age <- 0:maxage
  
  ## Age stuff
  nage <- length(age)
  msel <- rep(1,nage)
  
  ## for later use
  recruitmat <- matrix(0, nspace) # 4 is seasonality 
  recruitmat[1] <- 1 # 10 percent change of spawning north
  recruitmat[2] <- 1 # 90 percent of spawning south
  
  
  # Maturity
  
  movefifty <- movefiftyinit
  
  movemax <- rep(movemaxinit,nseason)
  
  movemat <- array(0, dim = c(nspace, nage, nseason, nyear)) # Chances of moving in to the other grid cell 
  
  
  if(nspace == 1){
    move = FALSE
  }else{
    move = TRUE
  }
  
  if(move == TRUE){
    for(j in 1:nspace){
      for(i in 1:nseason){
        movemat[j,,i,] <- movemax[i]/(1+exp(-moveslope*(age-movefifty)))
        
      }
    }
    
    movemat[,1:2,,] <- 0 # Recruits and 1 year olds don't move
    
    if(nseason == 4){ # For the standard model
      
      movemat[1,3:nage,2:3,] <- movesouth # Don't move south during the year
      movemat[1,3:nage,1,] <- moveout*0.5 # continuing south movement at spawning time
      
      movemat[1,3:nage,nseason,] <- moveout
      movemat[2,3:nage,nseason,] <- movesouth
    }
    # movemat[1,11:nage,nseason] <- 0
    # movemat[2,11:nage,nseason] <- 0
    
    
    
    # move.init <- array(0.5, dim = c(nspace, nage))
    # 
    # move.init[1,] <- 0.3
    # move.init[2,] <- 0.7
    move.init <- moveinit
    
  }else{
    move.init <- 1
  }
  
  # Run growth 
  Laa <- Linf*(1-exp(-K*(age-t0)))
  
  mat <- 1/(1+exp(-gamma*(age-tau)))
  wage_ssb <- wage_catch <- wage_survey <- wage_mid <- matrix(NA, nyear, nage)
  
  # Assign to matrices 
  
  wage_ssb <- t(replicate(nyear,0.01*Laa^3 ))
  wage_catch <- t(replicate(nyear,0.01*Laa^3 ))
  wage_survey <- t(replicate(nyear,0.01*Laa^3 ))
  wage_mid <- t(replicate(nyear,0.01*Laa^3 ))
 
  sel <- 1/(1+exp(-gamma_sel*(age-tau_sel)))
  
  if(length(catch) == 1){
    if(is.na(catch)){
      catch <- rep(NA, nyear)  
    }
  }
  
  if(fishing.type == 'constant'){
    if(length(F0) == 1){
      F0 <- rep(F0, nyear)
    }
    
  }
 
  
  if(fishing.type == 'AR'){
    
    Fin <- rep(0, nyear)
    Fin[1] <- F0
    omega <- rep(NA, nyear)
    omega[1] <- 0
    
    for(i in 2:nyear){
      
      omega[i] <- rho_F*omega[i-1]+sqrt(1-rho_F^2)*rnorm(1,0,sd = sqrt(SDF))
      Fin[i]<- Fin[1]*exp(omega[i]-0.5*(SDF^2))
      
      
    }
    F0 <- Fin
  }
  
  if(length(F0) != nyear){
    stop('wrong number of fishing mortality years')
  }
  
  # Option for natural mortality 
  if(mortality == 'constant' & length(M) == 1){
    M0 <-    rep(M, nyear)
  }
  
  if(mortality == 'rw'){ # Natural mortality random  walk (this needs a bias adjuster)
    M0 <- rep(0, nyear)
    M0[1] <- M
    
    for(i in 2:nyear){
    M0[i] <- M0[i-1]*exp(rnorm(n = 1, mean = 0, sd = SDM))
    
    if(M0[i] > M_limit){ # Specify a max limit for natural mortliaty 
      M0[i] <- M_limit
    }
    
    }
    
    
  }
  
  if(mortality == 'AR'){ # Natural mortality autocorrelation 
    
    M0 <- rep(0, nyear)
    M0[1] <- M
    omega <- rep(NA, nyear)
    omega[1] <- 0
    
    for(i in 2:nyear){
        
        omega[i] <- rho*omega[i-1]+sqrt(1-rho^2)*rnorm(1,0,sd = sqrt(SDM))
        M0[i]<-M0[1]*exp(omega[i]-0.5*(SDM^2))
        
        
    }
    
  }
  
  
  
  
  if(mortality == 'increase'){
    M0 <- rep(0, nyear)
    M0[1] <- M
    
    for(i in 2:nyear){
      M0[i] <- M0[i-1]+SDM
      
      if(M0[i] > M_limit){ # Specify a max limit for natural mortliaty 
        M0[i] <- M_limit
      }
      
    }
    
  }
  
  
  if(mortality == 'decrease'){
    
    M0 <- rep(0, nyear)
    M0[1] <- M
    
    for(i in 2:nyear){
      M0[i] <- M0[i-1]-SDM
      
      if(M0[i] < M_limit){ # Specify a max limit for natural mortliaty 
        M0[i] <- M_limit
      }
      M0[M0 < 0] <- 0
      
    }
  }
  
  
  if(length(M0) != nyear){
    stop('wrong number of natural mortality years')
  }

  
  if(nseason == 4){
  surveyseason <- 3
  
  }else{
    surveyseason <- floor(nseason/2)
  }
  
  if(nseason == 1){
    surveyseason = 1
  }
  
  # Insert Rdev function here 
  Rdev <- rnorm(nyear, mean = 0, sd = SDR)
    
  if(length(b) == 1){
  bout <- rep(b, nyear)
  }
  
  initN <- rep(0, nage-1) # 
  
  
  if(recruitment == 'BH'){
    R0 <- alpha/beta
  }
  
  parms <- list( # Keep parameters for an age based model in a seperate list
       logRinit = log(R0),
       logh = log(h),
       logSDsurv =log(SDsurv),
       #logSDR = log(1.4),
       #logphi_survey = log(11.33),
       # logSDF = log(0.1),
       # Selectivity parameters 
       initN = initN,
       Rin = Rdev     
       )
     
  flag_sel <- rep(0,nyear)
     
  if(nseason == 1){
    
    Fnseason <- matrix(rep(1, nspace))
    
  }else{
    
    if(exists('Fnseason') !=0)
    warning('Fisheries per season not provided. Assuming equal')
    Fnseason <- matrix(rep(1, nspace, nseason))/nseason
    
  }
  
 
  if(exists('omega') == FALSE){
    omega = NA
  }
  
  
  if(length(negg) == 1){
    negg <- rep(negg, nage)
  }
  
  if(length(negg) != nage){
    stop('wrong size')
  }
  
  egg.size <- (negg*wage_ssb[1,]^eggbeta)/1e6 # Scale to get easier numbers to work iwth
  
  
  
  df <-list(      #### Parameters #####
                  wage_ssb = t(wage_ssb),
                  wage_catch = t(wage_catch),
                  wage_survey = t(wage_survey),
                  wage_mid = t(wage_mid),
                  Msel = msel,
                  Matsel= as.numeric(mat),
                  nage = nage,
                  age = age,
                  nseason = nseason,
                  nyear = nyear,
                  tEnd = tEnd, 
                  logq = log(q),   
                  # Selectivity 
                  flag_sel = rep(0,nyear),
                  surveyseason = surveyseason,
                  M0 = M0,
                  K = K,
                  Linf = Linf,
                  nsurvey = nsurvey, # Frequency of survey years (e.g., 2 is every second year)
                  flag_survey = rep(1,nyear),
                  age_maxage = max(age), # Max age for age comps 
                  flag_catch =rep(1, nyear),
                  # variance parameters
                  logSDcatch = log(SDcatch),
                  logSDR = log(SDR), # Fixed in stock assessment ,
                  years = years,
                  b = bout,
                  bfuture = bfuture,
                  smul = smul, # Annual survey timing 
                  nspace = nspace,
                  Fpast = Fpast,
                  #TAC = TAC,
                  movemat = movemat,
                  move = move,
                  recruitmat = recruitmat,
                  move.init = move.init,
                  movefifty = movefifty,
                  movemax = movemax,
                  movesouth = movesouth,
                  moveout = moveout,
                  moveslope = moveslope,
                 # F0 = Fin,
                  sel = sel,
                  parms = parms,
                  Fnseason = Fnseason,
                  selectivity_change = selectivity_change,
                  Catch = catch,
                  catch.true = catch.true,
                  F0 = F0,
                  SDM = SDM,
                  rho = rho,
                  omega = omega,
                  M_limit = M_limit,
                  mortality = mortality,
                  recruitment = recruitment,
                  alpha = alpha,
                  beta = beta,
                  egg.size = egg.size
                  
                
                  # Parameters from the estimation model 
            )
  
  
  # Fix some things 
  if(df$age_maxage>max(df$age)){
    df$age_maxage <- max(df$age)
  }
  
  
  return(df)
  
}
