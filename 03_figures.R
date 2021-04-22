#load packages:
library("ggplot2")
library("cowplot")
library("reshape2")
library("dplyr")
library("gplots")

#set theme for plots
theme_set(theme_cowplot())


#load files:
tab <- read.csv2("\\\\path_to_file/global_data.csv") #raw data table
IC_tab <- read.csv2("\\\\path_to_file/IC_values_all_groups.csv") #produced in script 02

grandom_tab <- read.csv2("\\\\path_to_file/log_geo_fst_regression_randomized_and_original.csv") #produced in script 02
grandom_only=subset(grandom_tab,group=="random")
goriginal=subset(grandom_tab,group!="random")

erandom_tab <- read.csv2("\\\\path_to_file/e_div_fst_regression_randomized_and_original.csv") #produced in script 02
erandom_only=subset(erandom_tab,group=="random")
eoriginal=subset(erandom_tab,group!="random")

#calculate species / range median fst values:
tab$species_type <- paste(tab$species, tab$type, sep="_")
tab <- tab %>% group_by(species_type) %>% mutate(median_fst=median(fst))
tab.med <- tab[,c("type","parental_care","introduction_vector",
                  "single_multiple_introductions","feeding","median_fst",
                  "dispersal","species_type","gene","species","time_since_introduction_log")] %>% distinct()
tab.med$type.num <- as.numeric(tab.med$type=="non-native")+1 #native=1, non-native=2

############ PLOTS #################################################

#Figure 2:
tab.med$updown <- 0; tab.med$updown[tab.med$species %in% c("Podarcis siculus","Metrioptera roeselii")] <- 1

ggplot(tab,aes(x=type,y=fst,fill=type))+
  #geom_jitter(alpha=0.2,aes(col=type))+
  geom_boxplot(width=0.5)+
  geom_line(data=tab.med,aes(x=type.num,y=median_fst,group=species,col=as.factor(updown)),size=1.2,linetype="dotted")+
  geom_point(data=tab.med,aes(x=type,y=median_fst,col=as.factor(updown)),size=3.5)+
  scale_fill_manual(values=c("#5aae61","#762a83"))+
  scale_color_manual(values=c("black","darkgrey"),guide=FALSE)
#Figure was modified in inkscape.



#Figure 3A:
ggplot(tab,aes(x=introduction_vector,y=fst,fill=type))+
  geom_boxplot()+theme(legend.position='none')+
  scale_fill_manual(values=c("#5aae61","#762a83"))


#Figure 3B:
IC_tab2 <- melt(IC_tab,id.vars=c("group","N","chisq.p.fdr"),measure.vars=c("IC_ediv","IC_geo"))
IC_tab2$group[IC_tab2$group=="All"]="0All_0All"
IC_tab2$type <- as.factor(sapply(strsplit(as.character(IC_tab2$group),"\\_"), "[[", 1))
IC_tab2$type2 <- as.factor(sapply(strsplit(as.character(IC_tab2$group),"\\_"), "[[", 2))

ggplot(IC_tab2,aes(x=type2,y=value,fill=type, shape=variable))+
  geom_point(size=6,col="black")+
  scale_shape_manual(values=c(21,22))+
  scale_fill_manual(values=c("darkgrey","#5aae61","#762a83"))+
  coord_flip()
#Figure was modified in inkscape.






#Supplementary Figure 1:
ggplot(erandom_only,aes(x=t.value))+
  geom_histogram(binwidth=0.2)+
  geom_vline(data=eoriginal, aes(xintercept=t.value,col=group))

ggplot(grandom_only,aes(x=t.value))+
  geom_histogram(binwidth=0.2)+
  geom_vline(data=goriginal, aes(xintercept=t.value,col=group))
