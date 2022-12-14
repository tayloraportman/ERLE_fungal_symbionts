# Van Krevelin and Venn diagrams
#LCMS Isolates Alone

#Taylor Portman
#6NOV22

#Import libraries
library(tidyverse)
library(readxl)
library(ggpubr)
library(ggsci)
library(ggpolypath)
library(venn)
library(dplyr)
source('/Volumes/GoogleDrive/My Drive/projects/MastersThesis/R/ERLE_fungal_symbionts/LCMS/functions/funcitons_cdis_exploration.R')

#Set directories
project_dir <- setwd("/Volumes/GoogleDrive/My Drive/projects/MastersThesis/R/ERLE_fungal_symbionts/LCMS/IsolatesAlone")
figures_dir <- file.path(project_dir, paste0('RP_figures'))
tables_dir <- file.path(project_dir,  paste0('RP_data_clean'))
# Import data -------------------------------------------------------------
compounds_table<-read_csv('RP_data_clean/compounds_table.csv')
cd_results_table<-read.csv('RP_data_clean/gap_filled_compounds_table.csv')
no_gap_compounds_table<-read_csv('RP_data_clean/compounds_table_nogap.csv')

# Unique compounds --------------------------------------------------------
##Number of unique compounds (some Names/ compounds are replicated)
distinct<- distinct(cd_results_table, Name, .keep_all=TRUE)

# Van Krevelin -------------------------------------------------------------
# Ven Krevelen Diagram based on type of organic material
vank_material <- plot_vank(compounds_table, Class)
vank_material

figure_file <- file.path(figures_dir, 'Venk_diagram_all.jpg')
ggsave(figure_file, vank_material, width=30, height=20, unit= 'cm')

# Venn diagram ------------------------------------------------------------
# Venn diagram of the compounds present in each isolate
CTRL_list <- no_gap_compounds_table %>% 
  filter(Fungi == 'CTRL') %>% 
  select(FeatureID) %>% 
  distinct()

F1177_list <- no_gap_compounds_table %>% 
  filter(Fungi == '1177') %>% 
  select(FeatureID) %>% 
  distinct()

F1246_list <- no_gap_compounds_table %>% 
  filter(Fungi == '1246') %>% 
  select(FeatureID) %>% 
  distinct()

F1154_list <- no_gap_compounds_table %>% 
  filter(Fungi == '1154') %>% 
  select(FeatureID) %>% 
  distinct()

my_list <- list(CTRL = CTRL_list$FeatureID,
                F1177 = F1177_list$FeatureID,
                F1246 = F1246_list$FeatureID,
                F1154 = F1154_list$FeatureID)

my_colors <- c('CTRL' ='grey', 'F1154' = 'salmon', 'F1177' = 'coral4', 'F1246' = 'dodgerblue')

#Plot figure
venn_isolates <- plot_venn(my_list, my_colors)


#Save figure
figure_file <- file.path(figures_dir, 'venn_isolates2.png')
png(figure_file)
plot_venn <-  euler(my_list,
       quantities= TRUE,
       fills = my_colors,
       shape= 'ellipse',
       labels= list(font=4))
plot_venn
plot(plot_venn)
dev.off()





