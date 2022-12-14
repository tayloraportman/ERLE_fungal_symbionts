#Gina- Volcano plots from FTICR data

library(dplyr)
library(tidyverse)
#source('functions_cdis_diff.R')

data<- read.csv(file.choose())
dat<- data

metadataa<-read.csv (file.choose())
met<-metadataa

dat<- dat%>%
  mutate(FeatureID = paste0('Feature',formatC(n():0001, width = 4, flag = '0')))%>%
  gather(c(2:31), key = 'SampleID', value = 'AUC') %>% 
  filter(AUC > 0)

compounds_table <-  left_join(dat, met, by = 'SampleID')

########Get Samples######################################################################
get_samples <- function(metadata.df, Treatment, value){
  # Get value to filter samples
  selector <- syms({{Treatment}})
  
  samples <- metadata.df %>% 
    filter((!!! selector) == value)
  
  # Get only sample names
  samples <- samples$SampleID
  
  return(samples)
}
Drought.samples <- get_samples(met, Treatment = 'Treatment', value = 'Drought')
PreDrought.samples <- get_samples(met, Treatment = 'Treatment', value = 'PreDrought')
#F1246.samples <- get_samples(metadata, Treatment = 'Fungi', value = '1246')
#CTRL.samples <- get_samples(metadata, Treatment = 'Fungi', value = 'CTRL')
#####Diff table###############################################################################
get_diff_table <- function(auc_matrix, control.sample_list, treatment.sample_list, log2_transformed = FALSE){
  
  # Get the AUC values per each sample and calculate the means per feature
  temp.df_control <- auc_matrix %>% 
    select(all_of(control.sample_list))
  control_means <- rowMeans(temp.df_control, na.rm = TRUE)
  
  temp.df_treatment <- auc_matrix %>% 
    select(all_of(treatment.sample_list))
  treatment_means <- rowMeans(temp.df_treatment, na.rm = TRUE)
  
  diff_table <- as.data.frame(cbind(control_means, treatment_means))
  
  if(log2_transformed == TRUE){
    diff_table <- diff_table %>% 
      mutate(log2FC = treatment_means - control_means)
  } else {
    diff_table <- diff_table %>% 
      mutate(ratio = treatment_means/control_means) %>% # get the control/treatment ratio
      mutate(log2FC = log2(ratio)) # calculate log2FC
  }
  rownames(diff_table) <- rownames(auc_matrix)
  
  # Initialize pvalues matrix
  pvalues <- data.frame(row.names = rownames(auc_matrix), pval = rep(0, length(rownames(auc_matrix))))
  
  #Calculate pvalue per each of the features
  for(i in 1:nrow(pvalues)){
    t.test <- t.test(as.numeric(temp.df_control[i,]), as.numeric(temp.df_treatment[i,]), paired = FALSE)
    pvalues$pval[i] <- t.test$p.value
  }
  
  diff_table <- merge(diff_table, pvalues, by = 'row.names')
  diff_table$pval.adj <- p.adjust(diff_table$pval, method = 'fdr')
  diff_table <- diff_table %>%
    rename(FeatureID = Row.names)
  
  return(diff_table)
  
}

  
norm.matrix<- dat%>%
  select(FeatureID, SampleID, AUC)%>%
  pivot_wider('FeatureID',names_from='SampleID', values_from = 'AUC')
norm.matrix[is.na(norm.matrix)]<-0


Drought_to_PreDrought.diff_table <- get_diff_table(norm.matrix, treatment.sample_list = Drought.samples, control.sample_list = PreDrought.samples)      
#F1177_to_CTRL.diff_table <- get_diff_table(norm.matrix, treatment.sample_list = F1177.samples, control.sample_list = CTRL.samples)
#F1246_to_CTRL.diff_table <- get_diff_table(norm.matrix, treatment.sample_list = F1246.samples, control.sample_list = CTRL.samples)

# Create dataframes for the up and downregulated metabolites at each time point and merge them with the compound information
Drought_to_PreDrought.diff_table <- compounds_table %>% 
  select( - contains('Results'), -SampleID, -AUC,  -Treatment, - Plant) %>% 
  right_join(Drought_to_PreDrought.diff_table, by = 'FeatureID') %>% 
  distinct() %>% 
  filter(!is.nan(log2FC)) %>% 
  mutate(Comment = case_when(log2FC == Inf ~ 'Not present in control',
                             log2FC == -Inf ~ 'Only present in control',
                             log2FC < 0 ~ 'Downregulated',
                             log2FC > 0 ~ 'Upregulated'))

#F1154_to_CTRL.diff_table$Comment <- factor(F1154_to_CTRL.diff_table$Comment, 
                             #              levels = c('Only present in control', 'Downregulated', 'Upregulated', 'Not present in control' ))

#table_file <- file.path(tables_dir, 'Diff_expressed_F1154.csv')
#write_csv(F1154_to_CTRL.diff_table, table_file )


# Extract significant features of each comparison based on adjusted pvalue
sig_features <- c(Drought_to_PreDrought.diff_table$FeatureID[Drought_to_PreDrought.diff_table$pval < 0.05])
                #  T2_to_T0.diff_table$FeatureID[T2_to_T0.diff_table$pval < 0.05],
                 # T3_to_T0.diff_table$FeatureID[T3_to_T0.diff_table$pval < 0.05])
sig_features <- unique(sig_features)



# 4. Plots of dysregulated features


################################################################################################

## 4.2 Volcano plots
plot_volcano <- function(df, log2FC, pval, log2FC.threshold, pval.threshold){
  
  #Generate label for the plot
  
  ggplot(df,
         aes(x = {{log2FC}},
             y = -log10({{pval}}))) +
    geom_point(color = ifelse(abs(df$log2FC) > {{log2FC.threshold}} & 
                                -log10(df$pval) > -log10({{pval.threshold}}), "#FF0000", "#000000")) +
    geom_vline(xintercept = c(-{{log2FC.threshold}}, {{log2FC.threshold}}),
               linetype = 'dotted',
               size = 1,
               color = 'blue') +
    geom_hline(yintercept = -log10({{pval.threshold}}),
               linetype = 'dotted',
               size = 1,
               color = 'blue') +
    theme_bw() +
    labs(title = 'Volcano plot', 
         x = expression("Log"[2]*" Fold Change"), 
         y = expression("-Log"[10]*" pvalue")) +
    theme(plot.title = element_text(hjust = 0.5, 
                                    face = 'bold', 
                                    size = 18),
          plot.subtitle = element_text(hjust = 0.5, 
                                       face = 'bold', 
                                       size = 15))
}



lfc.t <- 2
pval.t <- 0.05
D_to_PD_volcano <- plot_volcano(Drought_to_PreDrought.diff_table, log2FC, pval, lfc.t, pval.t) +
  labs(subtitle = 'Drought vs PreDrought')
D_to_PD_volcano


####################################################IDK about this lastt part lol!!!!!#################

