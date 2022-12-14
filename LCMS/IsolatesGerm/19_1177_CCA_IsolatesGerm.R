#CCA
#Isolate 1177
#Isolates Germ

#TaylorPortman
#22OCT22

#Import libraries
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
norm.matrix<-read_csv('RP_data_clean/normalized_untransformed_auc_table.csv')
metadata<-read.csv("RP_data_clean/metadata.csv")
compounds_table<- read.csv('RP_data_clean/compounds_table.csv')


norm.matrix <- column_to_rownames(norm.matrix, var = '...1')
# Filter for only isolate 1177 --------------------------------------------
F1177.samples <- get_samples(metadata, Treatment = 'Fungi', value = '1177')
CTRL.samples<- c('41plus_RP_pos_MS2', '40plus_RP_pos_MS2','2plus_RP_pos_MS2')
CTRL.matrix<-norm.matrix%>%
  select(CTRL.samples)%>%
  rownames_to_column("FeatureID")
norm.matrix<-norm.matrix%>%
  select(F1177.samples)%>%
  rownames_to_column("FeatureID")%>%
  left_join(CTRL.matrix, by= 'FeatureID')%>%
  column_to_rownames(var="FeatureID")

norm.matrix<- t(norm.matrix)

CTRLmetadata<-metadata%>%
  filter(Fungi== 'CTRL', Plant== 'CTRL')
metadata<-metadata%>%
  filter(Fungi== '1177')%>%
  full_join(CTRLmetadata, by=c('SampleID', 'Plant', 'Fungi'))


SampleID = metadata[,1]
Plant = metadata[,2]
Fungi= metadata[,3]

SampleID = as.factor(SampleID)
Plant = as.factor(Plant)
Fungi = as.factor(Fungi)


CCA<- cca (formula= norm.matrix ~ Plant)
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
  geom_point(data=centroids, aes(CCA1, CCA2, color= c('black')), size=1)+
  geom_segment(data = centroids, aes(x = 0, y = 0, xend = (CCA1*1),
                                     yend = (CCA2*1)), arrow = arrow(length = unit(1/2, "picas")),
               color = "black")+
  geom_text_repel(data=centroids, 
                  x = (centroids$CCA1*1), 
                  y = (centroids$CCA2*1),
                  label=centroids$variable, max.overlaps = 12)+
  theme_classic()+
  labs(title= "Isolate 1177", caption= "Plant P = 0.01")
plot1
#Save to figures folder
fileName = paste(figures_dir, '1177_CCA_Plot.png',sep = '/')
ggsave(fileName, plot1, dpi = 800,width= 18, height=12, units=c('cm'))

#element data
plot2<- ggplot()+
  geom_point(data=sites, aes(CCA1, CCA2, color=Fungi, shape=Plant), size=3)+
  scale_color_manual(values=colors)+
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
fileName = paste(figures_dir, '1177_CCA_FeatureArrows.png',sep = '/')
ggsave(fileName, plot2, dpi = 800,width= 18, height=12, units=c('cm'))  
