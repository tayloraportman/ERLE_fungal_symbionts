###################################
#### New Volcano plot function ####
###################################
# Author: Christian Ayala

# -------------------------------------------------------------------------
plot_volcano <- function(df, log2FC, pval, log2FC.threshold, pval.threshold){
  
  #Generate label for the plot
  
  significant_points <- df %>% 
    select(FeatureID, {{log2FC}}, {{pval}}) %>% 
    filter(abs({{log2FC}}) > log2FC.threshold,
           -log10({{pval}}) > -log10(0.05)) %>% 
    pull(FeatureID)
  
  plot <- df %>%
    mutate(color4plot = ifelse(FeatureID %in% significant_points, 'significant', 'non-significant')) %>% 
    ggplot(aes(x = {{log2FC}},
               y = -log10({{pval}}))) +
    geom_point(aes(color = color4plot)) +
    scale_color_manual(values = c("significant" = 'red', 'non-significant' = 'black')) +
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
                                    face = 'bold'),
          plot.subtitle = element_text(hjust = 0.5,
                                       face = 'bold'),
          legend.position = 'none')
  
  return(plot)
}