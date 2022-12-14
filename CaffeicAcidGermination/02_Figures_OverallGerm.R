#Caffeic Acid Germination Experiment

#Overall germination figures
#Taylor Portman
#14OCT22

#Read in libraries
library(ggplot2)
library(dplyr)

#Set working directory
setwd("/Volumes/GoogleDrive/My Drive/projects/MastersThesis/R/ERLE_fungal_symbionts")

#Read in clean data
germ.dat<-read.csv("data_clean/CaffeicAcidGermination_Clean.csv")

#Figure save file path
figures= '/Volumes/GoogleDrive/My Drive/projects/MastersThesis/R/ERLE_fungal_symbionts/CaffeicAcidGermination/figures'

#Figure generation
ggplot(germ.dat, aes(Seed, Germ, fill=Media))+
  geom_boxplot()

#Remove DICA and CTRL from plot by removing from data
germ.dat2<-germ.dat%>%
  filter(Seed != "CTRL", Seed != "DICA")

plot<-ggplot(germ.dat2, aes(Seed, Germ, fill=Media))+
  geom_boxplot()+
  geom_point(position=position_dodge(width=0.75))+
  theme_bw()

#Save figure to figure file
fileName = paste(figures, 'AverageGerm_bySpecies.png',sep = '/')
ggsave(fileName, plot, dpi=400)


