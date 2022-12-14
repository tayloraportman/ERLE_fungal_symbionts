## Germination and mass loss for isolate 1154

#Taylor Portman
#March 3, 2022

library(dplyr)
library(ggplot2)
library(tidyr)
library(ggpubr)

germ<- read.csv("/Volumes/GoogleDrive/My Drive/projects/ERLE_fungal_symbionts/R/data/germination.csv")

#Take mean of controls for each Plant

# Plot germ for each (sample- Plant mean) as "Proportion of seeds germinated relative to control"
germ<- germ%>%
  select(plant_sp, fungal_id, plant_type, contam, final_div_initial)

control<- germ%>%
  filter(fungal_id== "CTRL2"| fungal_id=="CNTRL")%>%
  group_by(plant_sp) %>% summarize(mean = mean(final_div_initial))

germ<- germ%>%
  select(fungal_id, plant_sp, plant_type, final_div_initial)%>%

  filter(fungal_id=='1154')%>%
  mutate(diff = final_div_initial - control$mean)

means <- germ %>% 
  group_by(plant_sp, fungal_id, .add=TRUE) %>%
  summarise(
    means=mean(diff))


means<- means%>%
  mutate(plant_type = case_when(
    endsWith(plant_sp, "ERLE") ~ "Isolate + non-native plant",
    endsWith(plant_sp, "ERIN") ~ "Isolate + native plant",
    endsWith(plant_sp, "BOCU") ~ "Isolate + native plant",
    endsWith(plant_sp, "LEDU") ~ "Isolate + native plant"
  ))
germ<- germ%>%
  mutate(plant_type = case_when(
    endsWith(plant_sp, "ERLE") ~ "Isolate + non-native plant",
    endsWith(plant_sp, "ERIN") ~ "Isolate + native plant",
    endsWith(plant_sp, "BOCU") ~ "Isolate + native plant",
    endsWith(plant_sp, "LEDU") ~ "Isolate + native plant"
  ))

sd <- germ %>% 
  group_by(plant_sp, fungal_id, .add=TRUE) %>%
  summarise(
    sd=sd(diff))

means<- means%>%
  merge(left_join(means,sd,by=c("plant_sp","fungal_id")))


compare_means(diff ~ plant_sp, ref.group= "ERLE",  data = germ, method= 't.test')



plot<- ggplot(means, aes(x= plant_sp, y=means))+
  geom_bar(data=means, aes(x=plant_sp, y= means, fill= plant_type), alpha=0.7, stat='identity' )+
  geom_errorbar(data=means, aes(ymin= means-sd, ymax= means+sd), width=0.4, colour="black", alpha=0.7, size=0.5)+
  geom_point(data= germ, aes(x=plant_sp, y=diff, color=plant_type))+
  geom_hline(aes(yintercept= 0), color='black',linetype= 'dotted', size=1, )+
  scale_fill_manual(values=c("black", "red"))+
  scale_color_manual(values=c("black", "red"))+
  theme_classic() +
  guides(color = 'none')+
  labs(
    title = "Seed germination responce to ERLE endophyte 1154",
    fill= "Sample type")+
  xlab(label =  "Plant species") +
  ylab(label =  "Proportion of seeds germinated relative to controls")


plot+theme(text = element_text(size = 20))  


ggsave('1154_germination.png', plot, dpi = 800,width= 18, height=12, units=c('cm'))




##### Mass Loss #####
mass<- read.csv("/Volumes/GoogleDrive/My Drive/projects/ERLE_fungal_symbionts/R/data/mass.loss.csv")
mass$loss<-mass$X.Loss.Scaled.by.Controls
mass$plant_sp<- mass$Grass
mass$fungal_id<- as.factor(mass$Isolate)

mass<- mass%>%
  select(loss, plant_sp, fungal_id, plant_type)%>%
  filter(fungal_id== '1154')%>%
  mutate(plant_type = case_when(
    endsWith(plant_type, "Non-native") ~ "Isolate + non-native plant",
    endsWith(plant_type, "Native") ~ "Isolate + native plant"
  ))

#remove outliers by z score
mass.z<- mass %>% 
  mutate(zscore = (loss - mean(loss))/sd(loss))%>%
  filter(!(zscore>2.5| zscore<(-2.5)))

mean<- mass.z%>%
  group_by(plant_sp)%>%
  mutate(mean= (mean(loss)))

compare_means(loss ~ plant_sp,  data = mass.z, ref.group = "ERLE", method= 't.test')
my_comparisons<- list(c("BOCU","ERLE"), c("DICA", "ERLE"), c("SELE","ERLE"))

plot<- ggplot(mass.z, aes(x= plant_sp, y=loss, color=plant_type))+
  geom_boxplot(aes(x= plant_sp, y=loss,fill=plant_type), color='black',alpha=0.7)+
  geom_point(color='black')+
  #geom_hline(data=mean, aes(yintercept= mean, color=plant_type), linetype= 'dotted', size=1,alpha=0.7)+
  scale_fill_manual(values=c("black", "red"))+
  scale_color_manual(values=c("black", "red"))+
  stat_compare_means(comparisons= my_comparisons, ref.group = "ERLE", method= 't.test',label = "p.signif")+
 # stat_compare_means(method = "anova", label.y = 1.3)+ 
  guides(color = 'none')+
  theme_classic()+
  labs(
    title = "Litter decomposition percent mass loss by ERLE endophyte 1154",
    fill= "Sample type")+
  xlab(label =  "Plant species") +
  
  ylab(label =  "Percent mass loss relative to control")

plot+ theme(text = element_text(size = 40))    

ggsave('1154_mass.loss.png', plot, dpi = 800, width= 18, height=12, units=c('cm'))
