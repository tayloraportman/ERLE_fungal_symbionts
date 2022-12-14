#L2FC--Heatmaps Sig. only

#Taylor Portman
#7DEC22

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
project_dir <- setwd("/Volumes/GoogleDrive/My Drive/projects/MastersThesis/R/ERLE_fungal_symbionts/LCMS/IsolatesAlone")
figures_dir <- file.path(project_dir, paste0('RP_figures'))
tables_dir <- file.path(project_dir,  paste0('RP_data_clean'))
# Import data -------------------------------------------------------------
L2FCdat<- read_csv("RP_data_clean/L2FC_SigIsolates_BiogeochemicalClass.csv")
sig_all_compounds_table<-read_csv('RP_data_clean/Sig.diff.expressed.allIsolates_L2FC_0.5.csv')
norm.matrix<-read_csv('RP_data_clean/normalized_untransformed_auc_table.csv')
metadata<-read.csv("RP_data_clean/metadata.csv")
compounds_table<- read.csv('RP_data_clean/compounds_table.csv')





# Organize and combine data -----------------------------------------------
norm.matrix<- rename(norm.matrix, "FeatureID"="...1")

L2FCdat<-L2FCdat%>%
  select(FeatureID, Name, BiogeochemicalClass)%>%
  left_join(norm.matrix, by="FeatureID")%>%
  filter(!is.na(Name)) %>% 
  distinct()

#table_file <- file.path(tables_dir, 'Sig.diff.expressed.Filtered_biogeochemicalClass.csv')
#write_csv(L2FCdat, table_file )

metadata<-metadata %>%
  select(-Fungi2)


# Heatmap -----------------------------------------------------------------
col_annot <- column_to_rownames(metadata, var = 'SampleID') 
mapcolor <- colorRampPalette(brewer.pal(11, 'RdYlBu'))(100)[100:1]

names<- compounds_table%>%
  select(FeatureID, name4plot)



named_compounds<-L2FCdat%>%
  select(FeatureID, BiogeochemicalClass)%>%
  left_join(names, by="FeatureID")%>%
  distinct()
  
row_annot<- column_to_rownames(named_compounds, var='FeatureID')
row_annot<-row_annot%>%
  select(-FeatureID)


named_compounds.matrix <- norm.matrix[norm.matrix$FeatureID %in% named_compounds$FeatureID,]

named_compounds<- column_to_rownames(named_compounds, var = 'FeatureID')
named_compounds.matrix<-column_to_rownames(named_compounds.matrix, var= 'FeatureID')
row.names(named_compounds.matrix) <- compounds_table$name4plot[match(row.names(named_compounds.matrix), compounds_table$FeatureID)]


figure_file <- file.path(figures_dir, 'SigOnly_HeatmapNames.pdf')

annot_colors <- list(
  Fungi = c(CTRL ='grey', '1154' = 'salmon', '1177' = 'coral4', '1246' = 'dodgerblue'))

pdf(figure_file)
pheatmap(named_compounds.matrix,
         cluster_row_slices = TRUE,
         culstering_distance_rows= 'correlation',
         clustering_distance_cols = 'correlation',
         scale = 'row',
         annotation_col = col_annot,
         annotation_row= row_annot,
         annotation_colors = annot_colors,
         color = mapcolor,
         cutree_cols = 4,
         cutree_rows = 5,
         fontsize_row = 1,
         show_rownames = FALSE,
         main = 'Sig. features (scaled AUC)'
)
dev.off()



