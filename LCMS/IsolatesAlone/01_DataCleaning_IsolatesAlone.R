#Data Cleaning
#LCMS-- IsolatesAlone

#TaylorPortman
#6NOV22

#Import libraries
library(tidyverse)
library(readxl)
library(ggpubr)
library(dplyr)
source('/Volumes/GoogleDrive/My Drive/projects/MastersThesis/R/ERLE_fungal_symbionts/LCMS/functions/functions_cdis_exploration.R')

# Set up directories ------------------------------------------------------
project_dir <- setwd("/Volumes/GoogleDrive/My Drive/projects/MastersThesis/R/ERLE_fungal_symbionts/LCMS/IsolatesAlone")

# Give a name for your project so the results are all grouped together in a directory with the same name
project_name <- 'RP_IsolatesAlone'
figures_dir <- file.path(project_dir, paste0('RP_figures'))
tables_dir <- file.path(project_dir,  paste0('RP_data_clean'))

# Create output directories
#dir.create(figures_dir, showWarnings = FALSE)
#dir.create(tables_dir, showWarnings = FALSE)


# Import Metadata ---------------------------------------------------------
# Import metadata and fix names
metadata_file <- file.path(project_dir, 'data_raw', 'InputFiles_adj.csv')
metadata <- read_csv(metadata_file)
# Select only the useful columns and fix column names
metadata <- metadata %>% 
  select(SampleID,  Fungi) 

table_file <- file.path(tables_dir, 'metadata.csv')
write_csv(metadata, table_file)

# Import Data -------------------------------------------------------------
cd_results_file <- file.path(project_dir, 'data_raw', 'CompoundsChecked_IsolatesAlone.csv')
cd_results_table <- read_csv(cd_results_file)

# Naming and filtering ----------------------------------------------------
#REMOVE MANUALLY FILTERED DATA
cd_results_table<-cd_results_table%>%
  filter(Checked== "TRUE")
cd_results_table$DeltaMass<- as.numeric(cd_results_table$"Annot. DeltaMass [ppm]")
cd_results_table$MW<- as.numeric(cd_results_table$"Calc. MW")

cd_results_table<- cd_results_table%>%
  mutate(FeatureID = paste0('Feature',formatC(n():0001, width = 4, flag = '0'))) %>%
  select(FeatureID, Name, Formula, MW, contains('AnnotSource:'), contains('Results'), contains('Pathways'), 
        contains('Area:'), contains('Gap Status:')) %>% 
  # Differentiate between features that share the same name using "peak#" at the end of the name
  group_by(Name) %>%
  add_count(Name) %>% 
  # Create variable with names for plotting (useful in following scripts)
  mutate(name4plot = ifelse(is.na(Name), FeatureID, ifelse(n == 1, Name, paste0(Name, '-peak', n():1)))) %>% 
  dplyr:: select(!n) %>% 
  ungroup()

# Data manipulation -------------------------------------------------------
# Split formula column into elemental counts
cd_results_table <- separate_formula(cd_results_table)

# Calculate ratios and thermodynamic indices
cd_results_table <- calc_ratios_n_idxs(cd_results_table)

# Calculate classes
cd_results_table <- calc_classes(cd_results_table)

table_file <- file.path(tables_dir, 'RPcompounds.table.csv')
write_csv(cd_results_table, table_file)

# Gap filled Compounds Table ----------------------------------------------
# Gather area under the curve (AUC) values per sample
compounds_table <- cd_results_table %>% 
  select(-contains('Gap Status:')) %>% 
  gather(contains('Area:'), key = 'SampleID', value = 'AUC') %>% 
  filter(AUC > 0)
compounds_table$SampleID <- str_remove(compounds_table$SampleID, 'Area: ')

compounds_table$Class <- gsub( "Lignin", "Phenol", compounds_table$Class)
compounds_table$Class <- gsub("Tannin","Other",  compounds_table$Class)

# Save gap-filled table to be used in Statistical Analysis
table_file <- file.path(tables_dir, 'gap_filled_compounds_table.csv')
write_csv(compounds_table, table_file)

# Gap Filled Compounds Table ----------------------------------------------
# Gather area under the curve (AUC) values per sample
compounds_table <- cd_results_table %>% 
  select(-contains('Gap Status:')) %>% 
  gather(contains('Area:'), key = 'SampleID', value = 'AUC') %>% 
  filter(AUC > 0)
compounds_table$SampleID <- str_remove(compounds_table$SampleID, 'Area: ')

compounds_table$Class <- gsub( "Lignin", "Phenol", compounds_table$Class)
compounds_table$Class <- gsub("Tannin","Other",  compounds_table$Class)


# Save gap-filled table to be used in Statistical Analysis

table_file <- file.path(tables_dir, 'gap_filled_compounds_table.csv')
write_csv(compounds_table, table_file)

# No Gap Compounds Table  --------------------------------------------------------
# Gather "Gap Status" for filtering 
gap_status <- cd_results_table %>% 
  select(FeatureID, contains('Gap Status:')) %>% 
  gather(contains('Gap Status:'), key = 'SampleID', value = 'gap_status')
gap_status$SampleID <- str_remove(gap_status$SampleID, 'Gap Status: ')

## Filtering
compounds_table_nogap <- left_join(compounds_table, gap_status, by = c('FeatureID', 'SampleID')) 
compounds_table_nogap<-compounds_table_nogap%>%
  filter(gap_status!= "Full gap")

compounds_table <- left_join(compounds_table, gap_status, by = c('FeatureID', 'SampleID')) 

# Add metadata information
compounds_table_nogap$SampleID<- as.factor(compounds_table_nogap$SampleID)
metadata$SampleID<- as.factor(metadata$SampleID)
compounds_table_nogap <-  left_join(compounds_table_nogap, metadata, by = 'SampleID')


table_file <- file.path(tables_dir, 'compounds_table_nogap.csv')
write_csv(compounds_table_nogap, table_file)


# Compounds Table All --------------------------------------------------------------
# Add metadata information

compounds_table$SampleID<- as.factor(compounds_table$SampleID)
metadata$SampleID<- as.factor(metadata$SampleID)
compounds_table <-  left_join(compounds_table, metadata, by = 'SampleID')
compounds_table<- compounds_table


table_file <- file.path(tables_dir, 'compounds_table.csv')
write_csv(compounds_table, table_file)

