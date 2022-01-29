## Find the center of gravity of XX North Sea species ## 
library(tidyverse)
library(DATRAS)
library(surveyIndex)
#source('addSpatialData.R')

X <- readExchangeDir('IBTS/')
addSpatialData(X, "spatial data/ICES_areas.shp")
addSpectrum(X)
df <- read.csv('TEST.csv')

# Add latin names to Mikaels DF 

df$sciname <- NA
df$sciname[df$fishstock == "cod.27.47d20"] <- 'Gadus morhua'
df$sciname[df$fishstock == "had.27.46a20"] <- 'Melanogrammus aeglefinus'
df$sciname[df$fishstock == "her.27.3a47d"] <- "Clupea harengus"
df$sciname[df$fishstock == "nop.27.3a4"] <- 'Trisopterus esmarkii'
df$sciname[df$fishstock == "ple.27.420"] <- 'Pleuronectes platessa'
df$sciname[df$fishstock == "pok.27.3a46"] <- "Pollachius virens"
df$sciname[df$fishstock == "sol.27.4"] <- 'Solea solea'
df$sciname[df$fishstock == "spr.27.4"] <- 'Sprattus sprattus'
df$sciname[df$fishstock == "tur.27.4"] <- 'Scophthalmus maximus'
df$sciname[df$fishstock == "whg.27.47d"] <- 'Merlangius merlangus'
df$sciname[df$fishstock == "wit.27.3a47d"] <- 'Glyptocephalus cynoglossus'
df$sciname[df$fishstock == "spr.27.4"] <- 'Sprattus sprattus'
df$sciname[df$fishstock == "bss.27.4bc7ad-h"] <- 'Dicentrarchus labrax'
df$sciname[df$fishstock == "hke.27.3a46-8abd"] <- 'Merluccius merluccius'
df$sciname[df$fishstock == "hom.27.2a4a5b6a7a-ce-k8"] <- 'Trachurus trachurus'
df$sciname[df$fishstock ==  "mac.27.nea"] <- "Scomber scombrus"
df$sciname[df$fishstock == "san.4"] <- 'Ammodytes marinus' # Maybe add tobianus to this??

spp <- unique(df$sciname)

# Sort the data set by the assessed species # 

X <- subset.DATRASraw(X, Species %in% spp)

names(X$CA)
names(X$HH)
names(X$HL)

# Calculate the biomass in the haul 
X$CA$biomass <- X$CA$cp

# Create a data frame that combines haul id 

her <- subset.DATRASraw(X, Species == 'Clupea harengus')

str(her)

# Load the csv files 

csvs <- dir('IBTS_csv/')

for(i in 1:length(csvs)){
  
  xx <- read.csv(file.path('IBTS_csv',csvs[i]))
  
  if(i == 1){
    df.tot <- xx}
  else{
    df.tot <- rbind(df.tot, xx)
  }
  
}
df.tot$LngtClass <- df.tot$LngtClass*.1 # Convert to cm

df.tot$biomass <- .01*df.tot$LngtClass^3*df.tot$CPUE_number_per_hour

df.sum <- df.tot %>% group_by(Year, Species, LngtClass, ShootLat, ShootLong) %>%  
  summarise(biomass = mean(biomass)) %>% 
  group_by(Year, Species, ShootLat, ShootLong) %>% 
  summarise(cpue = sum(biomass))


# Calculate the biomass weighted mean distribution over time 
# THis is the simplest way to do it. Improve with VAST! 

df.cog <- df.sum %>% group_by(Year, Species) %>% 
  summarise(cogLat = weighted.mean(ShootLat, cpue),
            cogLong = weighted.mean(ShootLong, cpue))

p1 <- ggplot(df.cog, aes(x = Year, y = cogLat))+geom_line()+facet_wrap(~Species, scales = 'free_y')+
  scale_y_continuous('Latitude')+theme_bw()+geom_smooth()

p2 <- ggplot(df.cog, aes(x = Year, y = cogLong))+geom_line()+facet_wrap(~Species, scales = 'free_y')+
  scale_y_continuous('Latitude')+theme_bw()

print(p1)


print(p2)
