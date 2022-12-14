library (ggplot2)
library(mdthemes)

# Extract sample coordinates for PC1 and PC2
pca_coordinates <- as.data.frame(pca$x)
pca_coordinates$SampleID <- rownames(pca$x)

metadata<- metadata%>%
  mutate(PlantType = case_when(
    endsWith(Plant, "ERLE") ~ "Isolate + non-native plant",
    endsWith(Plant, "ERIN") ~ "Isolate + native plant",
    endsWith(Plant, "BOCU") ~ "Isolate + native plant",
    endsWith(Plant, "LEDU") ~ "Isolate + native plant",
    endsWith(Plant, "CTRL") ~ "Isolate alone"
  ))%>%
  mutate(Plant = case_when(
    endsWith(Plant, "ERLE") ~ "*E. lehmanniana* (ERLE)",
    endsWith(Plant, "ERIN") ~ "*E. intermedia* (ERIN)",
    endsWith(Plant, "BOCU") ~ "*B. curtipendula* (BOCU)",
    endsWith(Plant, "LEDU") ~ "*L. dubia* (LEDU)",
    endsWith(Plant, "CTRL") ~ "No plant (fungi alone)"
  ))%>%
  mutate(Fungi_name = case_when(
    endsWith(Fungi, "1154") ~ "*Fusarium sp.*(1154)",
    endsWith(Fungi, "1177") ~ "*Pseudophialophora sp.*(1177)",
    endsWith(Fungi, "1246") ~ "*Pseudothielavia sp.*(1246)",
    endsWith(Fungi, "CTRL") ~ "Media alone"
))
# Merge with metadata
pca <- left_join(pca_coordinates, metadata, by ='SampleID', na.rm=TRUE)
# Plot Individuals PCA
################################## PCA ALL##############################################
pca.plot <- ggplot(pca, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = Fungi_name,
                 shape = PlantType), size=4, alpha=0.7)+
  scale_color_manual( values=c('blue', 'darkgreen', 'darkorange', 'red'))  +
  scale_shape_manual( values=c(15,7,2,8))+
  theme_classic() +
  theme(plot.title = element_text(face = 'bold', hjust = 0.5))+
  mdthemes::md_theme_classic()+
  labs(
    title = "PCA ",
    color= "Sample type", shape= "Plant type", 
    caption= 'LC-MS/MS Isolate 1154')+
  xlab(label =  'PC1 (38.3%)') +
  ylab(label =  'PC2 (21.4%)')

pca.plot+theme(text = element_text(size = 20))  


figure_file <- file.path(figures_dir, 'PCA-revised2.png')

ggsave(figure_file, pca.plot, dpi = 800, width= 18, height=12, units=c('cm'))

##################################PCA ISOLATE 1154############################################
pca.plot <- ggplot(pca, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = PlantType,
                 shape = Plant), size=4, alpha=0.7)+
  scale_color_manual( values=c('black', 'red', 'darkgrey', 'grey'))  +
  scale_shape_manual( values=c(15,2,17,7,8))+
  theme_classic() +
  theme(plot.title = element_text(face = 'bold', hjust = 0.5))+
  mdthemes::md_theme_classic()+
labs(
  title = "PCA Isolate 1154",
  color= "Sample type", shape= "Plant species", 
  caption= 'LC-MS/MS Isolate 1154')+
  xlab(label =  'PC1 (41%)') +
  ylab(label =  'PC2 (19.5%)')

pca.plot+theme(text = element_text(size = 20))  


figure_file <- file.path(figures_dir, 'PCA-revised.png')
figure_file<- file.path( figures_dir, 'PCA-revised.tiff')
ggsave(figure_file, pca.plot, dpi = 800, width= 18, height=12, units=c('cm'))


# Prepare axis labels for PCA
pc1 <- paste0('PC1 (', round(eigen$variance.percent[1], digits = 1), '%)')
pc2 <- paste0('PC2 (', round(eigen$variance.percent[2], digits = 1), '%)')
