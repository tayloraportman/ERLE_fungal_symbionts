#SRER Endophyte Germination Assay

#Data Cleaning
#Taylor Portman
#14OCT22

#Read in libraries
library(dplyr)

#Set working directory
setwd("/Volumes/GoogleDrive/My Drive/projects/MastersThesis/R/ERLE_fungal_symbionts")

#Read in data
germ<-read.csv("data_raw/germination.csv")

#Save data_clean file path
data_clean = '/Volumes/GoogleDrive/My Drive/projects/MastersThesis/R/ERLE_fungal_symbionts/data_clean'


#Take mean of controls for each Plant

# Plot germ for each (sample- Plant mean) as "Proportion of seeds germinated relative to control"
germ<- germ%>%
  filter(! fungal_id== "1213")%>%
  dplyr:: select(plant_sp, fungal_id, plant_type, contam, final_div_initial )

control<- germ%>%
  filter(fungal_id== "CTRL2"| fungal_id=="CNTRL")%>%
  group_by(plant_sp) %>% summarize(mean = mean(final_div_initial))

germ<- germ%>%
  dplyr::select(fungal_id, plant_sp, plant_type, final_div_initial)%>%
  filter(!(fungal_id== "CTRL2"| fungal_id=="CNTRL"))%>%
  mutate(diff = final_div_initial - control$mean)

#Add means to each plant_sp
means <- germ %>% 
  group_by(plant_sp, fungal_id, .add=TRUE) %>%
  summarise(
    means=mean(diff))

sd <- germ %>% 
  group_by(plant_sp, fungal_id, .add=TRUE) %>%
  summarise(
    sd=sd(diff))

means<- means%>%
  merge(left_join(means,sd,by=c("plant_sp","fungal_id")))

dat<-means%>%
  merge(left_join(germ, means, by=c("plant_sp","fungal_id")))

### Rename columns and rows
dat<- dat%>%
  mutate(plant_type = case_when(
    endsWith(plant_sp, "ERLE") ~ "Isolate + non-native plant",
    endsWith(plant_sp, "ERIN") ~ "Isolate + native plant",
    endsWith(plant_sp, "BOCU") ~ "Isolate + native plant",
    endsWith(plant_sp, "LEDU") ~ "Isolate + native plant"
  ))%>%
  mutate(plant_sp = case_when(
    endsWith(plant_sp, "ERLE") ~ "*E. lehmanniana*",
    endsWith(plant_sp, "ERIN") ~ "*E. intermedia*",
    endsWith(plant_sp, "BOCU") ~ "*B. curtipendula*",
    endsWith(plant_sp, "LEDU") ~ "*L. dubia*",
    endsWith(plant_sp, "CTRL") ~ "No plant"
  ))

germ<- germ%>%
  mutate(plant_type = case_when(
    endsWith(plant_sp, "ERLE") ~ "Isolate + non-native plant",
    endsWith(plant_sp, "ERIN") ~ "Isolate + native plant",
    endsWith(plant_sp, "BOCU") ~ "Isolate + native plant",
    endsWith(plant_sp, "LEDU") ~ "Isolate + native plant"
  )) %>%
  mutate(plant_sp = case_when(
    endsWith(plant_sp, "ERLE") ~ "*E. lehmanniana*",
    endsWith(plant_sp, "ERIN") ~ "*E. intermedia*",
    endsWith(plant_sp, "BOCU") ~ "*B. curtipendula*",
    endsWith(plant_sp, "LEDU") ~ "*L. dubia*",
    endsWith(plant_sp, "CTRL") ~ "No plant"
  ))

# Export to data_clean ----------------------------------------------------

fileName = paste(data_clean, 'SREREndophyteGermination_Clean.csv',sep = '/')
write.csv(dat, fileName )
