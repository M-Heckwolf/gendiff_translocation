#load packages:
library("ggplot2")
library("cowplot")
library("reshape2")
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


############ PLOTS #################################################

#Figure 1:
ggplot(tab,aes(x=type,y=fst,fill=type))+
  geom_boxplot()+
  scale_fill_manual(values=c("#5aae61","#762a83"))



#Figure 2:
X=ggplot(tab,aes(x=dispersal,y=fst,fill=type))+
  geom_boxplot()+
  scale_fill_manual(values=c("#5aae61","#762a83"))

B=ggplot(tab,aes(x=introduction_vector,y=fst,fill=type))+
  geom_boxplot()+theme(legend.position='none')+
  scale_fill_manual(values=c("#5aae61","#762a83"))

Leg <- get_legend(X)

A <- X + theme(legend.position='none')

plot_grid(A, B, Leg, labels = c('A', 'B'), label_size = 14, ncol=3, rel_widths=c(1.5, 1.5,0.4))




#Figure 3:
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



#Supplementary Figure 2:
venn(list(active=rownames(tab[tab$dispersal=="active",]),passive=rownames(tab[tab$dispersal=="passive",]),intentional=rownames(tab[tab$introduction_vector=="intentional",]),unintentional=rownames(tab[tab$introduction_vector=="unintentional",])))
