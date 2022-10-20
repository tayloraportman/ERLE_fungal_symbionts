## 

#Taylor Portman
#Feb 26, 2022

library(dplyr)
library(ggplot2)
library(ggsignif)
library(devtools)
library(ggpubr)

getwd()
setwd("/Volumes/GoogleDrive/My Drive/projects/ERLE_fungal_symbionts/R")

mass<- read.csv("/Volumes/GoogleDrive/My Drive/projects/ERLE_fungal_symbionts/R/data/mass.loss.csv")
mass<-mass%>%
  filter(Isolate != "1213")
mass$loss<-mass$X.Loss.Scaled.by.Controls
mass$plant_sp<- mass$Grass
mass$fungal_id<- as.factor(mass$Isolate)

mass<- mass%>%
  dplyr::select(loss, plant_sp, fungal_id, plant_type)

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

compare_means(loss ~ plant_sp,  data = mass.z, ref.group = "ERLE", method= 't.test')
my_comparisons<- list(c("BOCU","ERLE"), c("DICA", "ERLE"), c("SELE","ERLE"))

plot<- ggplot(mass.z, aes(x= fungal_id, y=loss, color=plant_type))+
  geom_boxplot(aes(x= fungal_id, y=loss,fill=plant_type), color='black',alpha=0.7)+
  geom_point(color='black')+
  facet_grid(~plant_sp)+
  geom_hline(data=mean, aes(yintercept= mean, color=plant_type), linetype= 'dotted', size=1,alpha=0.7)+
scale_fill_manual(values=c("black", "red"))+
  scale_color_manual(values=c("black", "red"))+
  guides(color = 'none')+
  theme_classic()+
  labs(
    title = "Litter decomposition percent mass loss relative to control",
    fill= "Sample type",
    caption ="- - Mean by plant species")+
  xlab(label =  "Fungal isolate") +
  
  ylab(label =  "Percent mass loss relative to control")
  
plot+ theme(text = element_text(size = 40))    

ggsave('mass.loss.png', plot, dpi = 800, width= 18, height=12, units=c('cm'))


+ stat_compare_means(aes(group=plant_sp), label = "p.signif", method="t.test", comparisons = my_comparisons)
  







