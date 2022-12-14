#Log2FoldChange
#Isolates Germ

#TaylorPortman
#22OCT22

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
norm.matrix<-read_csv('RP_data_clean/normalized_untransformed_auc_table.csv')
metadata<-read.csv("RP_data_clean/metadata.csv")
compounds_table<- read.csv('RP_data_clean/compounds_table.csv')

norm.matrix <- column_to_rownames(norm.matrix, var = '...1')
# Calculate ratios for log2fold -------------------------------------------
# Get the average AUC per each of the treatments
## Get samples per treatment

## Fungal isolate treatment
F1154.samples <- get_samples(metadata, Treatment = 'Fungi', value = '1154')
F1177.samples <- get_samples(metadata, Treatment = 'Fungi', value = '1177')
F1246.samples <- get_samples(metadata, Treatment = 'Fungi', value = '1246')
CTRL.samples <- get_samples(metadata, Treatment = 'Fungi', value = 'CTRL')

# The following function will calculate ratio, log2FC, p values and adjusted pvalues. If no replicates are available for EACH treatment
# please use the get_diff_table_no_pval() function

F1154_to_CTRL.diff_table <- get_diff_table(norm.matrix, treatment.sample_list = F1154.samples, control.sample_list = CTRL.samples)
F1177_to_CTRL.diff_table <- get_diff_table(norm.matrix, treatment.sample_list = F1177.samples, control.sample_list = CTRL.samples)
F1246_to_CTRL.diff_table <- get_diff_table(norm.matrix, treatment.sample_list = F1246.samples, control.sample_list = CTRL.samples)


# Create dataframes for the up and downregulated metabolites at each time point and merge them with the compound information

F1154_to_CTRL.diff_table <- compounds_table %>% 
  select( - contains('Results'), -SampleID, -AUC, -Fungi) %>% 
  right_join(F1154_to_CTRL.diff_table, by = 'FeatureID') %>% 
  distinct() %>% 
  filter(!is.nan(log2FC)) %>% 
  mutate(Comment = case_when(log2FC == Inf ~ 'Not present in control',
                             log2FC == -Inf ~ 'Only present in control',
                             log2FC < 0 ~ 'Downregulated',
                             log2FC > 0 ~ 'Upregulated'))

F1154_to_CTRL.diff_table$Comment <- factor(F1154_to_CTRL.diff_table$Comment, 
                                           levels = c('Only present in control', 'Downregulated', 'Upregulated', 'Not present in control' ))

table_file <- file.path(tables_dir, 'Diff_expressed_F1154.csv')
write_csv(F1154_to_CTRL.diff_table, table_file )

F1177_to_CTRL.diff_table<- compounds_table %>% 
  select( - contains('Results'), -SampleID, -AUC, - Fungi) %>%
  right_join(F1177_to_CTRL.diff_table, by = 'FeatureID') %>% 
  distinct() %>% 
  filter(!is.nan(log2FC)) %>% 
  mutate(Comment = case_when(log2FC == Inf ~ 'Not present in control',
                             log2FC == -Inf ~ 'Only present in control',
                             log2FC < 0 ~ 'Downregulated',
                             log2FC > 0 ~ 'Upregulated'))

F1177_to_CTRL.diff_table$Comment <- factor(F1177_to_CTRL.diff_table$Comment, 
                                           levels = c('Only present in control', 'Downregulated', 'Upregulated', 'Not present in control' ))

table_file <- file.path(tables_dir, 'Diff_expressed_F1177.csv')
write_csv(F1177_to_CTRL.diff_table, table_file )

F1246_to_CTRL.diff_table <- compounds_table %>% 
  select(- contains('Results'), -SampleID, -AUC, -Fungi) %>% 
  right_join(F1246_to_CTRL.diff_table, by = 'FeatureID') %>% 
  distinct() %>% 
  filter(!is.nan(log2FC)) %>% 
  mutate(Comment = case_when(log2FC == Inf ~ 'Not present in control',
                             log2FC == -Inf ~ 'Only present in control',
                             log2FC < 0 ~ 'Downregulated',
                             log2FC > 0 ~ 'Upregulated'))

F1246_to_CTRL.diff_table$Comment <- factor(F1246_to_CTRL.diff_table$Comment, 
                                           levels = c('Only present in control', 'Downregulated', 'Upregulated', 'Not present in control' ))

table_file <- file.path(tables_dir, 'Diff_expressed_F1246.csv')
write_csv(F1246_to_CTRL.diff_table, table_file )


