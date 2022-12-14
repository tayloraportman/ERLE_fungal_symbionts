#RDA
#Normalized matrix
#LCMS Isolates Alone

#Taylor Portman
#10NOV22

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
source('/Volumes/GoogleDrive/My Drive/projects/MastersThesis/R/ERLE_fungal_symbionts/LCMS/functions/functions_cdis_norm_stats.R')

#Set directories
project_dir <- setwd("/Volumes/GoogleDrive/My Drive/projects/MastersThesis/R/ERLE_fungal_symbionts/LCMS/IsolatesAlone")
figures_dir <- file.path(project_dir, paste0('RP_figures'))
tables_dir <- file.path(project_dir,  paste0('RP_data_clean'))
# Import data -------------------------------------------------------------
norm.matrix<-read.csv('RP_data_clean/normalized_transformed_auc_table.csv', row.names=1)
metadata<-read.csv("RP_data_clean/metadata.csv")

SampleID = metadata[,1]
Fungi = metadata[,2]

SampleID = as.factor(SampleID)
Fungi = as.factor(Fungi)


CCA<- cca (formula= norm.matrix ~ Fungi)
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
  geom_point(data=sites, aes(CCA1, CCA2, color=Fungi), size=3)+
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

#element data
plot2<- ggplot()+
  geom_point(data=sites, aes(CCA1, CCA2, color=Fungi), size=3)+
 # scale_color_manual(values=colors)+
  geom_segment(data = elements, aes(x = 0, y = 0, xend = (CCA1*20),
                                    yend = (CCA2*20)), arrow = arrow(length = unit(1/2, "picas")),
               color = "black")+
  theme_bw() +
  geom_text_repel(data=elements, 
                  x = (elements$CCA1*20), 
                  y = (elements$CCA2*20),
                  label=elements$elements, max.overlaps = 12, size=3)+
  
  theme_classic()
plot2
#Save to figures folder
fileName = paste(figures, 'CCA_FeatureArrows.png',sep = '/')
ggsave(fileName, plot2, dpi = 800,width= 18, height=12, units=c('cm'))  
