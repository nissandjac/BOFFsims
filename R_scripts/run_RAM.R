# Run MSE from the RAM stocks # 
library(reshape2)
library(gtools)
library(patchwork)
library(tidyverse)
library(TMB)

source('R/calcSSB0.R')
source('R/run_agebased_model_true_Catch.R')
source('R/load_data_seasons.R')
source('R/load_data_future.R')

tsim <- 50
nyear <- 50
nruns <- 10
set.seed(1234)
saving <- TRUE
seeds <- round(runif(n = nruns, min = 1,  max = 1e6))

df <- load_data_seasons(nseason = 4, # Set up parameters 
                        Linf = 30, 
                        maxage = 10,
                        K = 1, 
                        t0 = 0, SDR = .7,
                        fishing.type = 'AR',
                        mortality = 'AR') # Specify parameters 

xrun <- run.agebased.true.catch(df)


plot(xrun$SSB)
