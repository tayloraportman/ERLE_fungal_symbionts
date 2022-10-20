## Germination of native and non native grasses grown in fungal endophytes from ERLE

#Taylor Portman
#Feb 26, 2022

library(dplyr)
library(ggplot2)
library(tidyr)
library(ggpubr)
library(mdthemes)

getwd()
setwd("/Volumes/GoogleDrive/My Drive/projects/ERLE_fungal_symbionts/R")

germ<- read.csv("/Volumes/GoogleDrive/My Drive/projects/ERLE_fungal_symbionts/R/data/germination.csv")

#Take mean of controls for each Plant

# Plot germ for each (sample- Plant mean) as "Proportion of seeds germinated relative to control"
germ<- germ%>%
 filter(! fungal_id== "1213")%>%
  dplyr:: select(plant_sp, fungal_id, plant_type, contam, final_div_initial )

control<- germ%>%
  filter(fungal_id== "CTRL2"| fungal_id=="CNTRL")%>%
  group_by(plant_sp) %>% summarize(mean = mean(final_div_initial))

germ<- germ%>%
  dplyr::select(fungal_id, plant_sp, plant_type, final_div_initial)%>%
  filter(!(fungal_id== "CTRL2"| fungal_id=="CNTRL"))%>%
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
  ))%>%
  mutate(plant_sp = case_when(
    endsWith(plant_sp, "ERLE") ~ "*E. lehmanniana*",
    endsWith(plant_sp, "ERIN") ~ "*E. intermedia*",
    endsWith(plant_sp, "BOCU") ~ "*B. curtipendula*",
    endsWith(plant_sp, "LEDU") ~ "*L. dubia*",
    endsWith(plant_sp, "CTRL") ~ "No plant"
  ))

germ<- germ%>%
  mutate(plant_type = case_when(
    endsWith(plant_sp, "ERLE") ~ "Isolate + non-native plant",
    endsWith(plant_sp, "ERIN") ~ "Isolate + native plant",
    endsWith(plant_sp, "BOCU") ~ "Isolate + native plant",
    endsWith(plant_sp, "LEDU") ~ "Isolate + native plant"
  )) %>%
  mutate(plant_sp = case_when(
    endsWith(plant_sp, "ERLE") ~ "*E. lehmanniana*",
    endsWith(plant_sp, "ERIN") ~ "*E. intermedia*",
    endsWith(plant_sp, "BOCU") ~ "*B. curtipendula*",
    endsWith(plant_sp, "LEDU") ~ "*L. dubia*",
    endsWith(plant_sp, "CTRL") ~ "No plant"
  ))

sd <- germ %>% 
  group_by(plant_sp, fungal_id, .add=TRUE) %>%
  summarise(
    sd=sd(diff))
 
means<- means%>%
  merge(left_join(means,sd,by=c("plant_sp","fungal_id")))





#### Stats####################################3#
compare_means(diff ~ fungal_id,  data = germ, method= 't.test') 
compare_means(diff ~ plant_sp, ref.group= "ERLE",  data = germ, method= 't.test')


#The experiment had very little sensitivity which has nothing to do with the poster.   
#You would need an MANOVA and some post-hoc test to see if the treatments had a statistically significant effect.

manova(diff~plant_sp, data=germ)

one<- aov(diff~ fungal_id, data=germ)
summary(one)

two<- aov(diff~ fungal_id + plant_sp, data=germ)
summary(two)

three<-aov(diff~ fungal_id * plant_sp, data=germ)
summary(three)

library(AICcmodavg)

model.set <- list(one, two, three)
model.names <- c("one", "two", "three")

aictab(model.set, modnames = model.names)


par(mfrow=c(2,2))
plot(two)
par(mfrow=c(1,1))

tukey.two.way<-TukeyHSD(two)

tukey.two.way

tukey.plot.aov<-aov(diff ~ fungal_id:plant_sp, data=germ)
tukey.plot.test<-TukeyHSD(tukey.plot.aov)
plot(tukey.plot.test, las = 1)

sub1<- subset(germ, fungal_id=="1154")
compare_means(diff ~ plant_sp,  data = sub1, method= 't.test') 
compare_means(diff ~ plant_sp, ref.group= "*E. lehmanniana*",  data = sub1, method= 't.test')

sub1<- subset(germ, fungal_id=="1177")
compare_means(diff ~ plant_sp,  data = sub1, method= 't.test') 
compare_means(diff ~ plant_sp, ref.group= "*E. lehmanniana*",  data = sub1, method= 't.test')

sub1<- subset(germ, fungal_id=="1246")
compare_means(diff ~ plant_sp,  data = sub1, method= 't.test') 
compare_means(diff ~ plant_sp, ref.group= "*E. lehmanniana*",  data = sub1, method= 't.test')

#### Plot ######################################

plot<- ggplot(means, aes(x= fungal_id, y=means))+
  geom_bar(data=means, aes(x=fungal_id, y= means, fill= plant_type), alpha=0.7, stat='identity' )+
  geom_errorbar(data=means, aes(ymin= means-sd, ymax= means+sd), width=0.4, colour="black", alpha=0.7, size=0.5)+
  geom_point(data= germ, aes(x=fungal_id, y=diff, color=plant_type))+
  geom_hline(aes(yintercept= 0), color='black',linetype= 'dotted', size=1, )+
  facet_grid(~plant_sp)+
  scale_fill_manual(values=c("black", "red"))+
  scale_color_manual(values=c("black", "red"))+
  theme_classic() +
  mdthemes::md_theme_classic() +
  guides(color = 'none')+
  labs(
    title = "Seed germination responce to ERLE endophytes",
    fill= "Sample type",
    caption ='Bars represent the mean number of seeds germinated relative to controls.')+
  xlab(label =  "Fungal isolate") +
  ylab(label =  "Proportion of seeds germinated relative to controls")+
  theme(axis.title.y = ggtext::element_markdown())


plot+theme(text = element_text(size = 20))  


ggsave('germination_metabo.png', plot, dpi = 800,width= 18, height=12, units=c('cm'))
  