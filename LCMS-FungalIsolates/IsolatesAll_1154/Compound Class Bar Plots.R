
## Compound Class Stacked Bar Plot
#For Ghiwa
#June 28, 2022


# I added a column to differentiate fungi between all the differentially expressed compounds


F1154.sig$Fungi<- '1154'
F1177.sig$Fungi<- '1177'
F1246.sig$Fungi<- '1246'


sig.data.all<- rbind(F1154.sig, F1177.sig, F1246.sig)
#This combines all my data sets
#You could also load your full diff. expressed compound list

#compounds_table$Class <- gsub( "Lignin", "Phenol", compounds_table$Class)
sig.data.all$Class <- gsub("Tannin","Other",  sig.data.all$Class)

#I used this to change some of the names of the classes
#In my final figure I added a different column in excll and manually 
#added compounds based on chem spider/ google searches
sig.data.all<- sig.data.all %>%
  mutate(class = case_when(
    endsWith(Class, "Lignin") ~ "Phenol",
    endsWith(Class, "Tannin") ~ "Highly oxygenated compounds"
  ))%>%
  select(FeatureID, Class, Fungi, Comment)%>%
  count(Class, Fungi, Comment)


########### PLOT 1--STACK BY TOTAL ABUNDANCE################
stack.class<-ggplot(sig.data.all, aes(fill=Class, y=n, x=Fungi)) + 
  facet_wrap(~Comment)+
  geom_bar(position="fill", stat="identity")+
  scale_fill_brewer(palette = 'Set2')+
  theme_classic()+
  ylab('Percent abundance of compounds')+
  xlab('Fungal isolates')+
  labs(title= "Differentially expressed compounds by compound class",
       comment="Differetally expressed compounds pval.adj< 0.05 ")
stack.class
figure_file <- file.path(figures_dir, 'Class.Percent.abund.png')
ggsave(figure_file, stack.class, dpi=400)
########### PLOT 2--STACK BY % ABUNDANCE################
stack.class2<-ggplot(sig.data.all, aes(fill=Class, y=n, x=Fungi)) + 
  facet_wrap(~Comment)+
  geom_bar(position="stack", stat="identity")+
  scale_fill_brewer(palette = 'Set2')+
  ylab('Number of compounds')+
  xlab('Fungal isolates')+
  labs(title= "Differentially expressed compounds by compound class",
       caption="Differetally expressed compounds pval< 0.05 ")+
  theme_classic()
stack.class2
figure_file <- file.path(figures_dir, 'Class.Relative.abund.png')
ggsave(figure_file, stack.class2, dpi=400)