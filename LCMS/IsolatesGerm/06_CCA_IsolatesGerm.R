#RDA
#Normalized matrix
#LCMS Isolates Germ

#Taylor Portman
#22NOV22

#Import libraries
library(readxl)
library(ggpubr)
library(ggsci)
library(gridExtra)
library(vegan)
library(tibble)
library(vegan)
library(factoextra)
library(vsn)
library(tidyverse)
library(ggplot2)
library(ggrepel)
source('/Volumes/GoogleDrive/My Drive/projects/MastersThesis/R/ERLE_fungal_symbionts/LCMS/functions/functions_cdis_norm_stats.R')

#Set directories
project_dir <- setwd("/Volumes/GoogleDrive/My Drive/projects/MastersThesis/R/ERLE_fungal_symbionts/LCMS/IsolatesGerm")
figures_dir <- file.path(project_dir, paste0('RP_figures'))
tables_dir <- file.path(project_dir,  paste0('RP_data_clean'))
# Import data -------------------------------------------------------------
norm.matrix<-read.csv('RP_data_clean/normalized_transformed_auc_table.csv', row.names=1)
metadata<-read.csv("RP_data_clean/metadata.csv")
compounds<- read.csv("RP_data_clean/RPcompounds.table.csv")

metadata2<-read.csv('data_raw/metadata.csv')

metadata2<-metadata2%>%
  select(SampleID, Plant, Fungi, Fungi2)



SampleID = metadata[,1]
Plant = metadata[,2]
Fungi= metadata[,3]
Fungi2= metadata2[,4]

SampleID = as.factor(SampleID)
Plant = as.factor(Plant)
Fungi = as.factor(Fungi)
Fungi2 = as.factor(Fungi2)


CCA<- cca (formula= norm.matrix ~ Fungi2+Plant)
CCA_stat<- anova.cca(CCA, by= "terms")

eigenvals(CCA)

plot_CCA<- plot(CCA, scaling=2)



# Ploting -----------------------------------------------------------------
elements<-as.data.frame(plot_CCA$species)
sites<-as.data.frame(plot_CCA$sites)
centroids<-as.data.frame(plot_CCA$centroids)

#add metadata 
sites<-sites%>%
  rownames_to_column(var="SampleID")%>%
  left_join(metadata, by="SampleID")
centroids<-centroids%>%
  rownames_to_column(var="variable")
elements<-elements%>%
  rownames_to_column(var="elements")


#PlotCentroids
colors <- c(CTRL ='grey', '1154' = 'salmon', '1177' = 'coral4', '1246' = 'dodgerblue')

plot1<-ggplot()+
  geom_point(data=sites, aes(CCA1, CCA2, color=Fungi, shape=Plant), size=3)+
scale_color_manual(values=colors)+
  geom_point(data=centroids, aes(CCA1, CCA2), size=1)+
  geom_segment(data = centroids, aes(x = 0, y = 0, xend = (CCA1*1),
                                     yend = (CCA2*1)), arrow = arrow(length = unit(1/2, "picas")),
               color = "black")+
  geom_text_repel(data=centroids, 
                  x = (centroids$CCA1*1), 
                  y = (centroids$CCA2*1),
                  label=centroids$variable, max.overlaps = 12)+
  theme_classic()
plot1
#Save to figures folder
fileName = paste(figures_dir, 'CCA_Plot.png',sep = '/')
ggsave(fileName, plot1, dpi = 800,width= 18, height=12, units=c('cm'))


# Calculate length of feature centroid from zero --------------------------

features<- elements%>%
  mutate(FeatureID= elements)%>%
  mutate(length= sqrt((CCA1*CCA1) + (CCA2*CCA2)))%>%
  filter(length> 0.07)

compounds<- compounds%>%
  select(FeatureID, Name)
  
features<- left_join(features, compounds, by= 'FeatureID')

#element data
plot2<- ggplot()+
  geom_point(data=sites, aes(CCA1, CCA2, color=Fungi, shape=Plant), size=3)+
  scale_color_manual(values=colors)+
  geom_point(data=features, aes(CCA1, CCA2), size=1)+
  geom_segment(data = features, aes(x = 0, y = 0, xend = (CCA1*10),
                                    yend = (CCA2*10)), arrow = arrow(length = unit(1/2, "picas")),
               color = "black")+
  theme_bw() +
  geom_text_repel(data=features, 
                  x = (features$CCA1*15), 
                  y = (features$CCA2*15),
                  label=features$Name, max.overlaps = 18, size=3)+
  
  theme_classic()
plot2
#Save to figures folder
fileName = paste(figures_dir, 'CCA_FeatureArrows.png',sep = '/')
ggsave(fileName, plot2, dpi = 800,width= 18, height=12, units=c('cm'))  
