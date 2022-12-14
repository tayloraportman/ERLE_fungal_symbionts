#UpsetPlots
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
library(tibble)
source('/Volumes/GoogleDrive/My Drive/projects/MastersThesis/R/ERLE_fungal_symbionts/LCMS/functions/functions_cdis_diff.R')

#Set directories
project_dir <- setwd("/Volumes/GoogleDrive/My Drive/projects/MastersThesis/R/ERLE_fungal_symbionts/LCMS/IsolatesGerm")
figures_dir <- file.path(project_dir, paste0('RP_figures'))
tables_dir <- file.path(project_dir,  paste0('RP_data_clean'))
# Import data -------------------------------------------------------------
ERLE_to_CTRL.diff_table<-read_csv('RP_data_clean/F1154_Diff_expressed_ERLE.csv')
ERIN_to_CTRL.diff_table<-read_csv('RP_data_clean/F1154_Diff_expressed_ERIN.csv')
BOCU_to_CTRL.diff_table<-read_csv('RP_data_clean/F1154_Diff_expressed_BOCU.csv')
LEDU_to_CTRL.diff_table<-read_csv('RP_data_clean/F1154_Diff_expressed_LEDU.csv')


# Extract sig features ----------------------------------------------------
# Extract significant features of each comparison based on adjusted pvalue
sig_features <- c(ERLE_to_CTRL.diff_table$FeatureID[ERLE_to_CTRL.diff_table$pval.adj < 0.05],
                  ERIN_to_CTRL.diff_table$FeatureID[ERIN_to_CTRL.diff_table$pval.adj < 0.05],
                  BOCU_to_CTRL.diff_table$FeatureID[BOCU_to_CTRL.diff_table$pval.adj < 0.05],
                  LEDU_to_CTRL.diff_table$FeatureID[LEDU_to_CTRL.diff_table$pval.adj < 0.05])

sig_features <- unique(sig_features)

# Number of different Features --------------------------------------------
num_diff_features <- tibble(Comparison = rep(c('ERLE_to_CTRL', 'ERIN_to_CTRL', 'BOCU_to_CTRL', 'LEDU_to_CTRL'), each = 2),
                            Type = rep(c('Name', 'No name'), 4),
                            count = 0)
num_diff_features$Type <- factor(num_diff_features$Type, levels = c('No name', 'Name'))

num_diff_features[1, 3] <- sum(ERLE_to_CTRL.diff_table$pval.adj < 0.05 & !is.na(ERLE_to_CTRL.diff_table$Name), na.rm = TRUE)
num_diff_features[2, 3] <- sum(ERLE_to_CTRL.diff_table$pval.adj < 0.05 & is.na(ERLE_to_CTRL.diff_table$Name), na.rm = TRUE)
num_diff_features[3, 3] <- sum(ERIN_to_CTRL.diff_table$pval.adj < 0.05 & !is.na(ERIN_to_CTRL.diff_table$Name), na.rm = TRUE)
num_diff_features[4, 3] <- sum(ERIN_to_CTRL.diff_table$pval.adj < 0.05 & is.na(ERIN_to_CTRL.diff_table$Name), na.rm = TRUE)
num_diff_features[5, 3] <- sum(BOCU_to_CTRL.diff_table$pval.adj < 0.05 & !is.na(BOCU_to_CTRL.diff_table$Name), na.rm = TRUE)
num_diff_features[6, 3] <- sum(BOCU_to_CTRL.diff_table$pval.adj < 0.05 & is.na(BOCU_to_CTRL.diff_table$Name), na.rm = TRUE)
num_diff_features[5, 3] <- sum(LEDU_to_CTRL.diff_table$pval.adj < 0.05 & !is.na(LEDU_to_CTRL.diff_table$Name), na.rm = TRUE)
num_diff_features[6, 3] <- sum(LEDU_to_CTRL.diff_table$pval.adj < 0.05 & is.na(LEDU_to_CTRL.diff_table$Name), na.rm = TRUE)


num_diff_plot <- ggplot(num_diff_features)+
  geom_bar(aes(y=count, x=Comparison,fill= Type), stat="identity", alpha=0.7)+
  theme_classic()+xlab("Comparison")+ylab("Count")

num_diff_plot
# shows snumber of significantly regulated compounds
figure_file <- file.path(figures_dir, 'Num_diff_features.png')
ggsave(figure_file, num_diff_plot, dpi = 300)


# Upset Plot ---------------------------------------------------------

#Filter diff tables by pval.adj

ERLE.sig<-ERLE_to_CTRL.diff_table%>% 
  filter(pval.adj<0.05)
ERIN.sig<-ERIN_to_CTRL.diff_table%>% 
  filter(pval.adj<0.05)
BOCU.sig<-BOCU_to_CTRL.diff_table%>% 
  filter(pval.adj<0.05)
LEDU.sig<-LEDU_to_CTRL.diff_table%>% 
  filter(pval.adj<0.05)

#Get list of features at each time point for each of the cases
ERLE.upregulated<-get_vectors(ERLE.sig, 'Comment', 'Upregulated', 'FeatureID')
ERLE.downregulated<-get_vectors(ERLE.sig, 'Comment', 'Downregulated', 'FeatureID')

ERIN.upregulated<-get_vectors(ERIN.sig, 'Comment', 'Upregulated', 'FeatureID')
ERIN.downregulated<-get_vectors(ERIN.sig, 'Comment', 'Downregulated', 'FeatureID')

BOCU.upregulated<-get_vectors(BOCU.sig, 'Comment', 'Upregulated', 'FeatureID')
BOCU.downregulated<-get_vectors(BOCU.sig, 'Comment', 'Downregulated', 'FeatureID')

LEDU.upregulated<-get_vectors(LEDU.sig, 'Comment', 'Upregulated', 'FeatureID')
LEDU.downregulated<-get_vectors(LEDU.sig, 'Comment', 'Downregulated', 'FeatureID')

upset_input<-list(ERLE.upregulated=ERLE.upregulated,
                  ERLE.downregulated=ERLE.downregulated,
                  
                  ERIN.upregulated=ERIN.upregulated,
                  ERIN.downregulated=ERIN.downregulated,
                  
                 BOCU.upregulated=BOCU.upregulated,
                  BOCU.downregulated=BOCU.downregulated,
                 
                 LEDU.upregulated=LEDU.upregulated,
                 LEDU.downregulated=LEDU.downregulated)

upset_metadata<-data.frame(sets=as.vector(names(upset_input)),
                           Time=rep(c('ERLE', 'ERIN', 'BOCU', 'LEDU'), each=4))

upset_features<-upset(fromList(upset_input), order.by="freq", cutoff=1, nsets=12,
                      mainbar.y.label='Number of shared features', sets.x.label='Number of features',
                      text.scale=c(1.5, 1, 1.5, 1, 1.3, 1.3),
                      set.metadata=list(data=upset_metadata, 
                                        plots=list(list(type='matrix_rows', column='Time', 
                                                        colors=c(ERLE='#EFC000FF', ERIN='#868686FF', BOCU='#CD534CFF', LEDU= '#003300')))))
upset_features

figure_file<-file.path(figures_dir, 'F1154_UpsetPlot.png')
png(figure_file, width=800, height=800, res=100)
upset_features
dev.off()