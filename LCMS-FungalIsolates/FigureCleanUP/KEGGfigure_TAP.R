#Roya Class composition up and down regulation figures

#F126.sig table is the significant up and down regulated 
#features which include comments for 
#Upregulation and Downregulation 

#Tables come from Chris' Differential Analysis script



## 4.3 KEGG Up and Down Regulated 


F1246.sig$KEGGPathways<-as.factor(F1246.sig$"KEGG Pathways")  

F1246.count<- F1246.sig%>% 
  select(FeatureID, KEGGPathways, Comment) %>% 
  separate_rows(KEGGPathways, sep = ';') %>% 
  group_by(KEGGPathways, Comment ) %>% 
  count(KEGGPathways) %>%
  drop_na(KEGGPathways)%>% filter(Comment != "Only present in control")
F1246.count <- within(F1246.count, n[Comment=="Downregulated"] <- (n*(-1)))
class.plot<- ggplot(F1246.count)+
  geom_bar(aes(y=KEGGPathways, x=n,fill= Comment), stat="identity", alpha=0.7)+
  scale_fill_manual(values = c("lightsteelblue","tomato2"))+
  theme_classic()+xlab("Number of Features")+ylab("KEGG Pathways")+
  labs(title="Up and down regulated KEGG pathways of Isolate 1246")
class.plot

figure_file <- file.path(figures_dir, 'KEGGrank_F1246.png')
ggsave(figure_file, class.plot, width=15, height=10, unit='cm')


F1177.sig$KEGGPathways<-as.factor(F1177.sig$"KEGG Pathways")  
F1177.count<- F1177.sig%>% 
  select(FeatureID, KEGGPathways, Comment) %>% 
  separate_rows(KEGGPathways, sep = ';') %>% 
  group_by(KEGGPathways, Comment ) %>% 
  count(KEGGPathways) %>%
  drop_na(KEGGPathways)%>% filter(Comment != "Only present in control")
F1177.count <- within(F1177.count, n[Comment=="Downregulated"] <- (n*(-1))) 
class.plot<- ggplot(F1177.count)+
  geom_bar(aes(y=KEGGPathways, x=n,fill= Comment), stat="identity", alpha=0.7)+
  scale_fill_manual(values = c("lightsteelblue","tomato2"))+
  theme_classic()+xlab("Number of Features")+ylab("KEGG Pathways")+
  labs(title="Up and down regulated KEGG pathways of Isolate 1177")
class.plot

figure_file <- file.path(figures_dir, 'KEGGrank_F1177.png')
ggsave(figure_file, class.plot, width=15, height=10, unit='cm')


F1154.sig$KEGGPathways<-as.factor(F1154.sig$"KEGG Pathways")  
F1154.count<- F1154.sig%>% 
  select(FeatureID, KEGGPathways, Comment) %>% 
  separate_rows(KEGGPathways, sep = ';') %>% 
  group_by(KEGGPathways, Comment ) %>% 
  count(KEGGPathways) %>%
  drop_na(KEGGPathways)%>% filter(Comment != "Only present in control")
F1154.count <- within(F1154.count, n[Comment=="Downregulated"] <- (n*(-1)))
class.plot<- ggplot(F1154.count)+
  geom_bar(aes(y=KEGGPathways, x=n,fill= Comment), stat="identity", alpha=0.7)+
  scale_fill_manual(values = c("lightsteelblue","tomato2"))+
  theme_classic()+xlab("Number of Features")+ylab("KEGG Pathways")+
  labs(title="Up and down regulated KEGG pathways of Isolate 1154")
class.plot

figure_file <- file.path(figures_dir, 'KEGGrank_F1154.png')
ggsave(figure_file, class.plot, width=15, height=10, unit='cm')


F1154.count$isolate<- 1154
F1177.count$isolate<- 1177
F1246.count$isolate<- 1246
class.rank.all<- rbind(F1154.count, F1177.count, F1246.count)

class.plot<- ggplot(class.rank.all)+
  geom_bar(aes(y=KEGGPathways, x=n,fill= Comment), stat="identity", alpha=0.7)+
  facet_wrap(~isolate)+
  scale_fill_manual(values = c("lightsteelblue","tomato2"))+
  theme_classic()+xlab("Number of Features")+ylab("KEGG Pathways")+
  labs(title="Up and down regulated KEGG pathways")
class.plot

figure_file <- file.path(figures_dir, 'KEGG_rank_all.png')
ggsave(figure_file, class.plot, width=30, height=10, unit='cm')

