# Data Normalization
#LC-MS/MS Isolates Alone

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
compounds_table<-read_csv('RP_data_clean/compounds_table.csv')
metadata<-read_csv('RP_data_clean/metadata.csv')


# Create AUC matrix -------------------------------------------------------
# Create a new tibble with the AUC per each mass from each sample
auc_table <- compounds_table %>% 
  select(FeatureID, SampleID, AUC)

# Transform the dataframe into a matrix-like table
auc_table <-  spread(auc_table, SampleID, AUC)
auc_table$FeatureID <- factor(auc_table$FeatureID, levels = str_sort(auc_table$FeatureID, numeric = TRUE))

auc_table <- auc_table %>% 
  arrange(FeatureID)
# Save untransformed data
auc_table <- column_to_rownames(auc_table, var = 'FeatureID')
table_file <- file.path(tables_dir, 'raw_auc_table.csv')
write.csv(auc_table, table_file, row.names = TRUE)


# Compare normalizaiton methods -------------------------------------------

normalization_plot <- normalize_by_all(auc_table)
figure_file <- file.path(figures_dir, 'all_normalized.boxplot.png')
ggsave(figure_file, normalization_plot, dpi = 300)


# Normalized matrix -------------------------------------------------------
# Obtaining transformed data for multivariate statistica analysis
norm.matrix <- cycloess.norm(auc_table)
table_file <- file.path(tables_dir, 'normalized_untransformed_auc_table.csv')
write.csv(norm.matrix, table_file, row.names= TRUE)

# Change missing values for zeroes
norm.matrix[is.na(norm.matrix)] <- 0

norm.matrix<-t(norm.matrix)
# Save normalized data
table_file <- file.path(tables_dir, 'normalized_transformed_auc_table.csv')
write.csv(norm.matrix, table_file, row.names= TRUE)

