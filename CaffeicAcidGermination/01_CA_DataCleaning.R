#Caffeic Acid Germination Experiment

#Data cleaning
#Taylor Portman
#14OCT22

#Read in libraries
library(tidyverse)

#Set workind directory
setwd("/Volumes/GoogleDrive/My Drive/projects/MastersThesis/R/ERLE_fungal_symbionts")

#Read in data
dat<-read.csv("data_raw/CaffeicAcid_SeedGermination.csv")

#Save data_clean file path
data_clean = '/Volumes/GoogleDrive/My Drive/projects/MastersThesis/R/ERLE_fungal_symbionts/data_clean'


###Data cleaning
#Rename day 14 data as total germination data
dat$Germ<-dat$Day14.30Sep
germ<-dat%>%
  dplyr::select( Sample, Seed, Media, Germ)

#add in means based on seed type (grass species)
means.grass <- germ %>% 
  dplyr::group_by(Seed, .add=TRUE) %>%
  dplyr::summarise(
    means=mean(Germ))

germ.dat<-germ%>%
  dplyr::left_join(means.grass, by="Seed", copy=FALSE)

#Normalize data using means
germ.dat$average<- germ.dat$Germ/ germ.dat$means

#Save data to data_clean
fileName = paste(data_clean, 'CaffeicAcidGermination_Clean.csv',sep = '/')
write.csv(germ.dat, fileName )
