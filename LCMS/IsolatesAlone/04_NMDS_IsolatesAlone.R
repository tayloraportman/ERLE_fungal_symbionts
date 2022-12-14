#NMDS
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
library(factoextra)
library(vsn)
library(tidyverse)
source('/Volumes/GoogleDrive/My Drive/projects/MastersThesis/R/ERLE_fungal_symbionts/LCMS/functions/functions_cdis_norm_stats.R')

#Set directories
project_dir <- setwd("/Volumes/GoogleDrive/My Drive/projects/MastersThesis/R/ERLE_fungal_symbionts/LCMS/IsolatesAlone")
figures_dir <- file.path(project_dir, paste0('RP_figures'))
tables_dir <- file.path(project_dir,  paste0('RP_data_clean'))
# Import data -------------------------------------------------------------
norm.matrix<-read.csv('RP_data_clean/normalized_transformed_auc_table.csv', row.names = 1)
metadata<-read.csv("RP_data_clean/metadata.csv")

#Set distance matrix method for relative abundance -ra calculations
nmds.matrix <- norm.matrix
dm.method <- 'bray'
dm <- vegdist(nmds.matrix, method=dm.method)

# NMDS  -------------------------------------------------------------------
color<- c( "#E64B35B2", "#4DBBD5B2",  "#00A087B2",  '#7E6148B2')
set.seed(123)
nmds <- metaMDS(dm,
                k = 2,
                maxit = 999,
                trymax = 500,
                wascores = TRUE)
stressplot(nmds)
# Extract nmds scores for plotting
nmds.scores <- as.data.frame(scores(nmds))
nmds.scores <- rownames_to_column(nmds.scores, var = 'SampleID')
nmds.scores <- left_join(nmds.scores, metadata, by = 'SampleID')
nmds_plot <- plot_nmds(nmds.scores, Fungi, color) +
  labs(title = 'NMDS plot by relative abundance')
nmds_plot
figure_file <- file.path(figures_dir, 'nmds_relative_abundance.png')
ggsave(figure_file, nmds_plot, dpi = 300)


# Permanova ---------------------------------------------------------------

rownames(metadata) <- metadata$SampleID
set.seed(456)
permanova <- adonis2(dm ~ Fungi, 
                     data= metadata, 
                     permutations=999, 
                     method="bray")
permanova #Marginally Significant (p= 0.086)

permanova <- adonis2(dm ~ Fungi2, #Fungi 2 groups isolates 1154 and 1177
                    data= metadata, 
                    permutations=999, 
                    method="bray")
permanova #Still marginally significant (p= 0.066)
