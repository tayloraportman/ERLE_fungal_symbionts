#SRER Endophyte Germination Assay

#Statistical Analysis
#Taylor Portman
#14OCT22

#Read in libraries
library(tidyverse)
library(AICcmodavg)
library(tidyr)
library(ggpubr)
library(mdthemes)

#Set working directory
setwd("/Volumes/GoogleDrive/My Drive/projects/MastersThesis/R/ERLE_fungal_symbionts")

#Read in data
germ<-read.csv("data_clean/SREREndophyteGermination_Clean.csv")

#Compare means with simple t.test
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

