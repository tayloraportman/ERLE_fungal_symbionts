##SRER Endophyte Germination Assay

#Plot figures
#Taylor Portman
#14OCT22

#Read in libraries
library(ggplot2)
library(ggpattern)
library(mdthemes)
library(dplyr)

#Set working directory
setwd("/Volumes/GoogleDrive/My Drive/projects/MastersThesis/R/ERLE_fungal_symbionts")

#Read in data
means<-read.csv("data_clean/SREREndophyteGermination_Clean.csv")
germ<-read.csv("data_clean/SREREndophyteGermination_Clean.csv")

means$fungal_id<-as.factor(means$fungal_id)
germ$fungal_id<- as.factor(germ$fungal_id)
means$plant_type<-as.factor(means$plant_type)
#Figure save file path
figures= 'ERLE_EndophyteGermination/figures'

#Plot Figures
colors <- c(CTRL ='grey', '1154' = 'salmon', '1177' = 'coral4', '1246' = 'dodgerblue')

plot<- ggplot(means, aes(x= fungal_id, y=means))+
  geom_bar(data=means, aes(x=fungal_id, y= means, fill= fungal_id), alpha=0.7, stat='identity' )+
  geom_errorbar(data=means, aes(ymin= means-sd, ymax= means+sd), width=0.4, color="black", size=0.5)+
  geom_point(data= germ, aes(x=fungal_id, y=diff, color=fungal_id))+
  geom_hline(aes(yintercept= 0), color='black', size=0.5)+
  facet_grid(~plant_sp)+
  scale_fill_manual(values=colors)+
  scale_color_manual(values=colors)+
  theme_classic() +
  mdthemes::md_theme_classic() +
  guides(color = 'none')+
  labs(
    title = "Seed germination responce to ERLE endophytes",
    fill= "Sample type",
    caption ='Bars represent the mean number of seeds germinated relative to controls.')+
  xlab(label =  "Fungal isolate") +
  ylab(label =  "Proportion of seeds germinated relative to controls")+
  theme(axis.title.y = ggtext::element_markdown())


plot+theme(text = element_text(size = 20))  

# Save to figures folder
fileName = paste(figures, 'EndophyteGerm_Final.png',sep = '/')
ggsave(fileName, plot, dpi = 800,width= 18, height=12, units=c('cm'))
