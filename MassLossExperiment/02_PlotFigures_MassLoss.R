#Litter Mass Loss In Vitro

#Plot figures
#Taylor Portman
#14OCT22

#Read in libraries
library(ggplot2)
library(ggsignif)
library(devtools)
library(ggpubr)

#Set working directory
setwd("/Volumes/GoogleDrive/My Drive/projects/MastersThesis/R/ERLE_fungal_symbionts")

#Read in data
mean<-read.csv("data_clean/InVitro_MassLoss_Clean.csv")

#Figure save file path
figures= '/Volumes/GoogleDrive/My Drive/projects/MastersThesis/R/ERLE_fungal_symbionts/MassLossExperiment/figures'

#comparisons
compare_means(loss ~ plant_sp,  data = mass.z, ref.group = "ERLE", method= 't.test')
my_comparisons<- list(c("BOCU","ERLE"), c("DICA", "ERLE"), c("SELE","ERLE"))

#Plot figure: mass loss
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
#theme(text = element_text(size = 40))   
#stat_compare_means(aes(group=plant_sp), label = "p.signif", method="t.test", comparisons = my_comparisons)

plot
fileName = paste(figures, 'InVitro_LitterMassLoss_Final.png',sep = '/')
ggsave(fileName, plot, dpi = 800,width= 18, height=12, units=c('cm'))
