#
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
sig_all_compounds_table<-read_csv('RP_data_clean/Sig.diff.expressed.allIsolates_L2FC_0.5.csv')


# VanKrevlin --------------------------------------------------------------
my_colors=c('Upregulated'='tomato2',
            'Downregulated'='lightsteelblue')
dysreg_vk<-ggplot(sig_all_compounds_table,
                  aes(x=O_to_C,
                      y=H_to_C,
                      color=Class)) +
  geom_point() +
  scale_color_manual(values=my_colors) +
  theme_classic() +
  labs(x='O:C',
       y='H:C',
       title='Van Krevelen Diagram ',
       color='Assigned compound class') +
  #new_scale_color() +
  #class_rect+
  scale_color_manual(values=get_palette('d3', 8)) +
  theme(plot.title=element_text(face='bold', hjust=0.5)) +
  facet_wrap(~Fungi*Comment, 3)

dysreg_vk


# Compound Class ----------------------------------------------------------
my_colors=c('Upregulated'='tomato2',
            'Downregulated'='lightsteelblue')

count.all<- sig_all_compounds_table%>% 
  select(FeatureID, Class, Comment, pval.adj, pval, Fungi)%>%
  count(Fungi,Class, Comment)

count.all <- within(count.all, n[Comment=="Downregulated"] <- (n*(-1)))

class.plot<- ggplot(count.all)+
  geom_bar(aes(y=Class, x= n, fill= Comment), stat="identity", alpha=0.7)+
  facet_wrap(~Fungi, 1)+
  scale_fill_manual(values = my_colors)+
  theme_classic()+xlab("Number of Features")+ylab("Class")+
  labs(title="Class ranks")
class.plot

figure_file <- file.path(figures_dir, 'CompoundClass_rank_all.png')
ggsave(figure_file, class.plot, width=15, height=10, unit='cm')

class.plot<- ggplot(count.all)+
  geom_bar(aes(y=Class, x=n, fill= Comment), stat="identity", alpha=0.7)+
  facet_wrap(~Fungi, 3)+
  scale_fill_manual(values = c("lightsteelblue","tomato2"))+
  theme_classic()+xlab("Number of Features")+ylab("Class")+
  labs(title="Class ranks")
class.plot

figure_file <- file.path(figures_dir, 'CompoundClass_rank_all_stacked.png')
ggsave(figure_file, class.plot, width=15, height=10, unit='cm')



# GFE of sig. compounds ---------------------------------------------------------------------

all_GFE_box<- ggplot(sig_all_compounds_table,
                     aes(x = Comment,
                         y = GFE,
                         fill = Comment)) +
  geom_violin() +
  geom_boxplot(width=0.1, color="black", alpha=0.2) +
  scale_fill_manual(values = my_colors) +
  theme_classic() +
  theme(plot.title = element_text(face = 'bold', hjust = 0.5))+
  facet_wrap(~Fungi)
all_GFE_box

figure_file<-file.path(figures_dir, 'all_GFE.png')
ggsave(figure_file, all_GFE_box)








