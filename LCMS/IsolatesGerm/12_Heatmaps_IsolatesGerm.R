#Heatmaps/ Hierarchical Clustering
#IsolatesGerm

#TaylorPortman
#10OCT22

#Import libraries
library(tidyverse)
library(ggrepel)
library(ggsci)
library(pheatmap)
library(RColorBrewer)
library(readxl)
library(UpSetR)
source('/Volumes/GoogleDrive/My Drive/projects/MastersThesis/R/ERLE_fungal_symbionts/LCMS/functions/functions_cdis_diff.R')

#Set directories
project_dir <- setwd("/Volumes/GoogleDrive/My Drive/projects/MastersThesis/R/ERLE_fungal_symbionts/LCMS/IsolatesGerm")
figures_dir <- file.path(project_dir, paste0('RP_figures'))
tables_dir <- file.path(project_dir,  paste0('RP_data_clean'))
# Import data -------------------------------------------------------------
sig_all_compounds_table<-read_csv('RP_data_clean/Sig.diff.expressed.allIsolates_L2FC_0.5.csv')
norm.matrix<-read_csv('RP_data_clean/normalized_untransformed_auc_table.csv')
metadata<-read.csv("RP_data_clean/metadata.csv")
compounds_table<- read.csv('RP_data_clean/compounds_table.csv')

norm.matrix <- column_to_rownames(norm.matrix, var = '...1')

# Heatmap of all detected features ----------------------------------------
# Initialize graphical device
# Set sampleID as row.names to annotate heatmap

col_annot <- column_to_rownames(metadata, var = 'SampleID') 
mapcolor <- colorRampPalette(brewer.pal(11, 'RdYlBu'))(100)[100:1]

figure_file <- file.path(figures_dir, 'All_features_heatmap.pdf')

annot_colors <- list(
  Fungi = c(CTRL ='#0073C2FF', F1154 = '#EFC000FF', F1177 = '#868686FF', F1246 = '#CD534CFF'),
  Plants= c(ERLE= 'goldenrod', ERIN= 'purple', BOCU= 'green', LEDU= 'blue'))

pdf(figure_file)
pheatmap(norm.matrix,
         clustering_distance_rows = 'correlation',
         clustering_distance_cols = 'correlation',
         scale = 'row',
         annotation_col = col_annot,
       #  annotation_colors = annot_colors,
         color = mapcolor,
         show_rownames = FALSE,
         cutree_cols = 5,
         main = 'All features (scaled AUC)'
)
dev.off()



# Extract named features --------------------------------------------------
# Extract all features that have names
named_compounds <- compounds_table %>% 
  select(FeatureID, Name, name4plot) %>% 
  filter(!is.na(Name)) %>% 
  distinct()

# Create a matrix of only identified features

named_compounds.matrix <- norm.matrix[rownames(norm.matrix) %in% named_compounds$FeatureID,]
row.names(named_compounds.matrix) <- compounds_table$name4plot[match(row.names(named_compounds.matrix), compounds_table$FeatureID)]

# ID features Heatmap -----------------------------------------------------
figure_file <- file.path(figures_dir, 'All_identified_features_heatmap.pdf')

pdf(figure_file, width = 15, height = 10)
pheatmap(named_compounds.matrix,
         clustering_distance_rows = 'correlation',
         clustering_distance_cols = 'correlation',
         scale = 'row',
         annotation_col = col_annot,
         #annotation_colors = annot_colors,
         color = mapcolor,
         cutree_cols = 4,
         #cutree_rows = 5,
         #fontsize_row = 6,
        show_rownames = FALSE,
         main = 'Identified features (scaled AUC)'
)
dev.off()
