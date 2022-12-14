#Extract Significant Features
#IsolatesAlone

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
project_dir <- setwd("/Volumes/GoogleDrive/My Drive/projects/MastersThesis/R/ERLE_fungal_symbionts/LCMS/IsolatesAlone")
figures_dir <- file.path(project_dir, paste0('RP_figures'))
tables_dir <- file.path(project_dir,  paste0('RP_data_clean'))
# Import data -------------------------------------------------------------
compounds_table<- read.csv('RP_data_clean/compounds_table.csv')
F1154_to_CTRL.diff_table<-read_csv('RP_data_clean/Diff_expressed_F1154.csv')
F1177_to_CTRL.diff_table<-read_csv('RP_data_clean/Diff_expressed_F1177.csv')
F1246_to_CTRL.diff_table<-read_csv('RP_data_clean/Diff_expressed_F1246.csv')

# Filter diff tables by pval.adj
F1246.sig <- F1246_to_CTRL.diff_table %>% 
  filter(pval.adj < 0.05)
F1154.sig <- F1154_to_CTRL.diff_table %>% 
  filter(pval.adj < 0.05)
F1177.sig <- F1177_to_CTRL.diff_table %>% 
  filter(pval.adj < 0.05)
# Get list of features at each time point for each of the cases
F1246.upregulated <- get_vectors(F1246.sig, 'Comment', 'Upregulated', 'FeatureID')
F1246.downregulated <- get_vectors(F1246.sig, 'Comment', 'Downregulated', 'FeatureID')
F1154.upregulated <- get_vectors(F1154.sig, 'Comment', 'Upregulated', 'FeatureID')
F1154.downregulated <- get_vectors(F1154.sig, 'Comment', 'Downregulated', 'FeatureID')
F1177.upregulated <- get_vectors(F1177.sig, 'Comment', 'Upregulated', 'FeatureID')
F1177.downregulated <- get_vectors(F1177.sig, 'Comment', 'Downregulated', 'FeatureID')



# Significant Features Table ----------------------------------------------

sig_all<-rbind(F1154.sig%>% select(FeatureID, Name, Comment, log2FC, pval, pval.adj) %>% mutate(Fungi='1154'),
               F1177.sig%>% select(FeatureID, Name, Comment, log2FC, pval, pval.adj) %>% mutate(Fungi='1177'),
               F1246.sig%>% select(FeatureID, Name, Comment, log2FC, pval, pval.adj) %>% mutate(Fungi='1246'))

sig_all_compounds_table<-compounds_table%>% 
  select(FeatureID, H_to_C, O_to_C, Class, Fungi, GFE) %>% 
  left_join(sig_all, by=c('FeatureID', 'Fungi'))%>%
  distinct() %>% 
  na.omit()

table_file <- file.path(tables_dir, 'Sig.diff.expressed.allIsolates_L2FC_0.5.csv')
write_csv(sig_all_compounds_table, table_file )

