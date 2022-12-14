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
source('/Volumes/GoogleDrive/My Drive/projects/MastersThesis/R/ERLE_fungal_symbionts/LCMS/functions/functions_cdis_diff.R')

#Set directories
project_dir <- setwd("/Volumes/GoogleDrive/My Drive/projects/MastersThesis/R/ERLE_fungal_symbionts/LCMS/IsolatesGerm")
figures_dir <- file.path(project_dir, paste0('RP_figures'))
tables_dir <- file.path(project_dir,  paste0('RP_data_clean'))
# Import data -------------------------------------------------------------
F1154_to_CTRL.diff_table<-read_csv('RP_data_clean/Diff_expressed_F1154.csv')
F1177_to_CTRL.diff_table<-read_csv('RP_data_clean/Diff_expressed_F1177.csv')
F1246_to_CTRL.diff_table<-read_csv('RP_data_clean/Diff_expressed_F1246.csv')



# Extract sig features ----------------------------------------------------
# Extract significant features of each comparison based on adjusted pvalue
sig_features <- c(F1154_to_CTRL.diff_table$FeatureID[F1154_to_CTRL.diff_table$pval.adj < 0.05],
                  F1177_to_CTRL.diff_table$FeatureID[F1177_to_CTRL.diff_table$pval.adj < 0.05],
                  F1246_to_CTRL.diff_table$FeatureID[F1246_to_CTRL.diff_table$pval.adj < 0.05])

sig_features <- unique(sig_features)

# Number of different Features --------------------------------------------
num_diff_features <- tibble(Comparison = rep(c('F1154_to_CTRL', 'F1177_to_CTRL', 'F1246_to_CTRL'), each = 2),
                            Type = rep(c('Name', 'No name'), 3),
                            count = 0)
num_diff_features$Type <- factor(num_diff_features$Type, levels = c('No name', 'Name'))

num_diff_features[1, 3] <- sum(F1154_to_CTRL.diff_table$pval.adj < 0.05 & !is.na(F1154_to_CTRL.diff_table$Name), na.rm = TRUE)
num_diff_features[2, 3] <- sum(F1154_to_CTRL.diff_table$pval.adj < 0.05 & is.na(F1154_to_CTRL.diff_table$Name), na.rm = TRUE)
num_diff_features[3, 3] <- sum(F1177_to_CTRL.diff_table$pval.adj < 0.05 & !is.na(F1177_to_CTRL.diff_table$Name), na.rm = TRUE)
num_diff_features[4, 3] <- sum(F1177_to_CTRL.diff_table$pval.adj < 0.05 & is.na(F1177_to_CTRL.diff_table$Name), na.rm = TRUE)
num_diff_features[5, 3] <- sum(F1246_to_CTRL.diff_table$pval.adj < 0.05 & !is.na(F1246_to_CTRL.diff_table$Name), na.rm = TRUE)
num_diff_features[6, 3] <- sum(F1246_to_CTRL.diff_table$pval.adj < 0.05 & is.na(F1246_to_CTRL.diff_table$Name), na.rm = TRUE)


num_diff_plot <- ggplot(num_diff_features)+
  geom_bar(aes(y=count, x=Comparison,fill= Type), stat="identity", alpha=0.7)+
  theme_classic()+xlab("Comparison")+ylab("Count")

num_diff_plot
# shows snumber of significantly regulated compounds
figure_file <- file.path(figures_dir, 'Num_diff_features.png')
ggsave(figure_file, num_diff_plot, dpi = 300)


# Upset Plot ---------------------------------------------------------

#Filter diff tables by pval.adj

F1154.sig<-F1154_to_CTRL.diff_table%>% 
  filter(pval.adj<0.05)
F1177.sig<-F1177_to_CTRL.diff_table%>% 
  filter(pval.adj<0.05)
F1246.sig<-F1246_to_CTRL.diff_table%>% 
  filter(pval.adj<0.05)

#Get list of features at each time point for each of the cases
F1154.upregulated<-get_vectors(F1154.sig, 'Comment', 'Upregulated', 'FeatureID')
F1154.downregulated<-get_vectors(F1154.sig, 'Comment', 'Downregulated', 'FeatureID')

F1177.upregulated<-get_vectors(F1177.sig, 'Comment', 'Upregulated', 'FeatureID')
F1177.downregulated<-get_vectors(F1177.sig, 'Comment', 'Downregulated', 'FeatureID')

F1246.upregulated<-get_vectors(F1246.sig, 'Comment', 'Upregulated', 'FeatureID')
F1246.downregulated<-get_vectors(F1246.sig, 'Comment', 'Downregulated', 'FeatureID')


upset_input<-list(F1154.upregulated=F1154.upregulated,
                  F1154.downregulated=F1154.downregulated,
                  
                  F1177.upregulated=F1177.upregulated,
                  F1177.downregulated=F1177.downregulated,
                
                  F1246.upregulated=F1246.upregulated,
                  F1246.downregulated=F1246.downregulated)

upset_metadata<-data.frame(sets=as.vector(names(upset_input)),
                           Time=rep(c('F1154', 'F1177', 'F1246'), each=4))

upset_features<-upset(fromList(upset_input), order.by="freq", cutoff=1, nsets=12,
                      mainbar.y.label='Number of shared features', sets.x.label='Number of features',
                      text.scale=c(1.5, 1, 1.5, 1, 1.3, 1.3),
                      set.metadata=list(data=upset_metadata, 
                                        plots=list(list(type='matrix_rows', column='Time', 
                                                        colors=c(F1154='#EFC000FF', F1177='#868686FF', F1246='#CD534CFF')))))
upset_features

figure_file<-file.path(figures_dir, 'UpsetPlot.png')
png(figure_file, width=800, height=800, res=100)
upset_features
dev.off()

