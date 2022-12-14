#PCA
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
library(factoextra)
library(vsn)
library(tidyverse)
source('/Volumes/GoogleDrive/My Drive/projects/MastersThesis/R/ERLE_fungal_symbionts/LCMS/functions/functions_cdis_norm_stats.R')

#Set directories
project_dir <- setwd("/Volumes/GoogleDrive/My Drive/projects/MastersThesis/R/ERLE_fungal_symbionts/LCMS/IsolatesGerm")
figures_dir <- file.path(project_dir, paste0('RP_figures'))
tables_dir <- file.path(project_dir,  paste0('RP_data_clean'))
# Import data -------------------------------------------------------------
norm.matrix<-read.csv('RP_data_clean/normalized_transformed_auc_table.csv', row.names= 1)
metadata<-read.csv("RP_data_clean/metadata.csv")


color<- c( "#E64B35B2", "#4DBBD5B2",  "#00A087B2",  '#7E6148B2')

# PCA ---------------------------------------------------------------------
# Calculate PCA with prcomp
pca <- prcomp(norm.matrix)
# Get eigenvalues
eigen <- get_eigenvalue(pca)
# Plot screeplot using the functions from factoextra
scree_plot <- fviz_eig(pca, addlabels = TRUE) +
  theme_bw() +
  theme(plot.title = element_text(face = 'bold', hjust = 0.5))
scree_plot 
figure_file <- file.path(figures_dir, 'screeplot.png')
ggsave(figure_file, scree_plot, dpi = 300)
# Plot cumulative variance plot
cumvar_plot <- plot_cumvar(eigen)
cumvar_plot
figure_file <- file.path(figures_dir, 'cumulative_variance.png')
ggsave(figure_file, cumvar_plot, dpi = 300)

# Plot PCA ----------------------------------------------------------------
# Extract sample coordinates for PC1 and PC2'
metadata$SampleID<- as.factor(metadata$SampleID)
pca_coordinates <- as.tibble(pca$x)
pca_coordinates$SampleID <- rownames(pca$x)
# Merge with metadata
pca_coordinates <- left_join(pca_coordinates, metadata, by ='SampleID')
# Prepare axis labels for PCA
pc1 <- paste0('PC1 (', round(eigen$variance.percent[1], digits = 1), '%)')
pc2 <- paste0('PC2 (', round(eigen$variance.percent[2], digits = 1), '%)')
# Plot Individuals PCA
pca_plot <- plot_dotplot(pca_coordinates, PC1, PC2, Fungi, Plant) +
  labs(title = 'PCA plot',
       x = pc1,
       y = pc2)
pca_plot
figure_file <- file.path(figures_dir, 'PCA-plot.png')
ggsave(figure_file, pca_plot, dpi = 300)

