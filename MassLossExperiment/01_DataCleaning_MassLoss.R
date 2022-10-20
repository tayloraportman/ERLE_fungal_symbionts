#Litter Mass Loss In Vitro

#Data Cleaning
#Taylor Portman
#14OCT22

#Read in libraries
library(dplyr)

#Set working directory
setwd("/Volumes/GoogleDrive/My Drive/projects/MastersThesis/R/ERLE_fungal_symbionts")

#Read in data
mass<-read.csv("data_raw/mass.loss.csv")

#Save data_clean file path
data_clean = '/Volumes/GoogleDrive/My Drive/projects/MastersThesis/R/ERLE_fungal_symbionts/data_clean'


#Data clean up
mass<-mass%>%
  filter(Isolate != "1213")
mass$loss<-mass$X.Loss.Scaled.by.Controls
mass$plant_sp<- mass$Grass
mass$fungal_id<- as.factor(mass$Isolate)

mass<- mass%>%
  select(loss, plant_sp, fungal_id, plant_type)

#remove outliers by z score
mass.z<- mass %>% 
  mutate(zscore = (loss - mean(loss))/sd(loss))%>%
  filter(!(zscore>2.5| zscore<(-2.5)))%>%
  mutate(plant_type = case_when(
    endsWith(plant_type, "Non-native") ~ "Isolate + non-native plant",
    endsWith(plant_type, "Native") ~ "Isolate + native plant"
  ))

mean<- mass.z%>%
  group_by(plant_sp)%>%
  mutate(mean= (mean(loss)))

fileName = paste(data_clean, 'InVitro_MassLoss_Clean.csv',sep = '/')
write.csv(mean, fileName )


