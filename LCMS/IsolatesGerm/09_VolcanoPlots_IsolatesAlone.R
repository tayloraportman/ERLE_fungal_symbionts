#Volcano Plots
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
F1154_to_CTRL.diff_table<-read_csv('RP_data_clean/Diff_expressed_F1154.csv')
F1177_to_CTRL.diff_table<-read_csv('RP_data_clean/Diff_expressed_F1177.csv')
F1246_to_CTRL.diff_table<-read_csv('RP_data_clean/Diff_expressed_F1246.csv')


# Volcano Plot -------------------------------------------------------------


lfc.t <- 0.3
pval.t <- 0.05

F1154_to_CTRL_volcano <- plot_volcano(F1154_to_CTRL.diff_table, log2FC, pval.adj, lfc.t, pval.t) +
  labs(subtitle = 'Isolate 1154 vs Media')
F1154_to_CTRL_volcano

figure_file <- file.path(figures_dir, 'F1154_CTRL_volcano.png')
ggsave(figure_file, F1154_to_CTRL_volcano, width=20, height=20, unit='cm')


F1177_to_CTRL_volcano <- plot_volcano(F1177_to_CTRL.diff_table, log2FC, pval.adj, lfc.t, pval.t) +
  labs(subtitle = 'Isolate 1177 vs Media')
F1177_to_CTRL_volcano

figure_file <- file.path(figures_dir, 'F1177_CTRL_volcano.png')
ggsave(figure_file, F1177_to_CTRL_volcano, width=20, height=20, unit='cm')

F1246_to_CTRL_volcano <- plot_volcano(F1246_to_CTRL.diff_table, log2FC, pval.adj, lfc.t, pval.t) +
  labs(subtitle = 'Isolate 1246 vs Media')
F1246_to_CTRL_volcano

figure_file <- file.path(figures_dir, 'F1246_CTRL_volcano.png')
ggsave(figure_file, F1246_to_CTRL_volcano, width=20, height=20, unit='cm')


