#Log2FoldChange
#Isolate 1177
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

CTRLmetadata<-metadata%>%
  filter(Fungi== 'CTRL', Plant== 'CTRL')
metadata<-metadata%>%
  filter(Fungi== '1177')%>%
  full_join(CTRLmetadata, by=c('SampleID', 'Plant', 'Fungi'))

#Select treatments= Plant 
ERLE.samples<-get_samples(metadata, Treatment= 'Plant', value= 'ERLE')
ERIN.samples<-get_samples(metadata, Treatment= 'Plant', value= 'ERIN')
BOCU.samples<-get_samples(metadata, Treatment= 'Plant', value= 'BOCU')
LEDU.samples<-get_samples(metadata, Treatment= 'Plant', value= 'LEDU')
CTRL.plant.samples<-get_samples(metadata, Treatment= 'Plant', value= 'CTRL')

ERLE_to_CTRL.diff_table<- get_diff_table(norm.matrix, treatment.sample_list = ERLE.samples, control.sample_list = CTRL.plant.samples)
ERIN_to_CTRL.diff_table<- get_diff_table(norm.matrix, treatment.sample_list = ERIN.samples, control.sample_list = CTRL.plant.samples)
BOCU_to_CTRL.diff_table<- get_diff_table(norm.matrix, treatment.sample_list = BOCU.samples, control.sample_list = CTRL.plant.samples)
LEDU_to_CTRL.diff_table<- get_diff_table(norm.matrix, treatment.sample_list = LEDU.samples, control.sample_list = CTRL.plant.samples)


# Create dataframes for the up and downregulated metabolites at each time point and merge them with the compound information

ERLE_to_CTRL.diff_table <- compounds_table %>% 
  filter(Plant== "ERLE")%>%
  select( - contains('Results'), -gap_status, -SampleID, -AUC, -Fungi, -Plant) %>% 
  right_join(ERLE_to_CTRL.diff_table, by = 'FeatureID') %>% 
  distinct() %>% 
  filter(!is.nan(log2FC)) %>% 
  mutate(Comment = case_when(log2FC == Inf ~ 'Not present in control',
                             log2FC == -Inf ~ 'Only present in control',
                             log2FC < 0 ~ 'Downregulated',
                             log2FC > 0 ~ 'Upregulated'))

ERLE_to_CTRL.diff_table$Comment <- factor(ERLE_to_CTRL.diff_table$Comment, 
                                          levels = c('Only present in control', 'Downregulated', 'Upregulated', 'Not present in control' ))

table_file <- file.path(tables_dir, 'F1177_Diff_expressed_ERLE.csv')
write_csv(ERLE_to_CTRL.diff_table, table_file )

#ERIN
ERIN_to_CTRL.diff_table <- compounds_table %>% 
  filter(Plant== "ERIN")%>%
  select( - contains('Results'), -gap_status, -SampleID, -AUC, -Fungi, -Plant) %>% 
  right_join(ERIN_to_CTRL.diff_table, by = 'FeatureID') %>% 
  distinct() %>% 
  filter(!is.nan(log2FC)) %>% 
  mutate(Comment = case_when(log2FC == Inf ~ 'Not present in control',
                             log2FC == -Inf ~ 'Only present in control',
                             log2FC < 0 ~ 'Downregulated',
                             log2FC > 0 ~ 'Upregulated'))

ERIN_to_CTRL.diff_table$Comment <- factor(ERIN_to_CTRL.diff_table$Comment, 
                                          levels = c('Only present in control', 'Downregulated', 'Upregulated', 'Not present in control' ))

table_file <- file.path(tables_dir, 'F1177_Diff_expressed_ERIN.csv')
write_csv(ERIN_to_CTRL.diff_table, table_file )

#BOCU
BOCU_to_CTRL.diff_table <- compounds_table %>% 
  filter(Plant== "BOCU")%>%
  select( - contains('Results'), -gap_status, -SampleID, -AUC, -Fungi, -Plant) %>% 
  right_join(BOCU_to_CTRL.diff_table, by = 'FeatureID') %>% 
  distinct() %>% 
  filter(!is.nan(log2FC)) %>% 
  mutate(Comment = case_when(log2FC == Inf ~ 'Not present in control',
                             log2FC == -Inf ~ 'Only present in control',
                             log2FC < 0 ~ 'Downregulated',
                             log2FC > 0 ~ 'Upregulated'))

BOCU_to_CTRL.diff_table$Comment <- factor(BOCU_to_CTRL.diff_table$Comment, 
                                          levels = c('Only present in control', 'Downregulated', 'Upregulated', 'Not present in control' ))

table_file <- file.path(tables_dir, 'F1177_Diff_expressed_BOCU.csv')
write_csv(BOCU_to_CTRL.diff_table, table_file )

#LEDU
LEDU_to_CTRL.diff_table <- compounds_table %>% 
  filter(Plant== "LEDU")%>%
  select( - contains('Results'), -gap_status, -SampleID, -AUC, -Fungi, -Plant) %>% 
  right_join(LEDU_to_CTRL.diff_table, by = 'FeatureID') %>% 
  distinct() %>% 
  filter(!is.nan(log2FC)) %>% 
  mutate(Comment = case_when(log2FC == Inf ~ 'Not present in control',
                             log2FC == -Inf ~ 'Only present in control',
                             log2FC < 0 ~ 'Downregulated',
                             log2FC > 0 ~ 'Upregulated'))

LEDU_to_CTRL.diff_table$Comment <- factor(LEDU_to_CTRL.diff_table$Comment, 
                                          levels = c('Only present in control', 'Downregulated', 'Upregulated', 'Not present in control' ))

table_file <- file.path(tables_dir, 'F1177_Diff_expressed_LEDU.csv')
write_csv(LEDU_to_CTRL.diff_table, table_file )
