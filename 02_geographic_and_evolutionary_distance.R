#load packages
library("ggplot2")
library("cowplot")
library("hier.part")

#set theme for plots
theme_set(theme_cowplot())

#load table
tab <- read.csv2("\\\\path_to_file/global_data.csv")


#------------------> geographic and evolutionary distances
#we were interested in the association between geographic / evolutionary distance and genetic differentiation. 
#Therefore, we fitted a linear regression using FST as the dependent and log-transformed geographic or evolutionary distances as the explanatory variable. 
#In the next step, the dataset was split up into eight sub-groups, which were chosen according to their significant differences in FST values 
#(interaction effect of type x introduction vector; Table 2). 
#These analysis will give us the absolute vaiation explained by evolutionary / geographic distance (R2 values)
#Additionally, in order to compare the relative proportion of variation in FST explained by evolutionary and geographic distance 
#we used the hierarchical partitioning approach (IC values).




############### absolute variation explained by evolutionary distance (e_div)
#we found two interesting interactions that we want to zoom into:
#type x introduction vector

#we want to know whether e_div explains more variation in fst 
#    in one group compared to the other within the interaction groups

#create the groups you want to look at:
tab$type_introvect <- paste(tab$type,tab$introduction_vector,sep="_")


#we are interested in the question: How much variation in fst is explained by e_div?
#Thus regression is a good model to use
#data is not normally distributed, thats why we cannot trust the p value calculated for the regression
#with a bootstrap randomization test we create our own distribution
#
#randomization tests (10,000x bootstrap)
efit_original <- lm(fst~e_div,data=tab)

#create table with original and randomized values:
erandom_tab=data.frame(group="All",
                       R2=signif(summary(efit_original)$r.squared, 5),
                       Adj_R2=signif(summary(efit_original)$adj.r.squared, 5),
                       Intercept=signif(efit_original$coef[[1]],5 ),
                       Slope=signif(efit_original$coef[[2]], 5),
                       Estimate=signif(summary(efit_original)$coef[2,1], 5),
                       Std.Error=signif(summary(efit_original)$coef[2,2], 5),
                       t.value=signif(summary(efit_original)$coef[2,3], 5),
                       P=signif(summary(efit_original)$coef[2,4], 5),
                       N=dim(tab)[1])

###add the values for the 4 groups (introduction*type)
for (i in unique(tab$type_introvect)){
  edata = subset(tab[tab$type_introvect==i,])
  efit1 <- lm(fst~e_div,data=edata)
  erandom_tab=rbind(erandom_tab,data.frame(group=i,
                                           R2=signif(summary(efit1)$r.squared, 5),
                                           Adj_R2=signif(summary(efit1)$adj.r.squared, 5),
                                           Intercept=signif(efit1$coef[[1]],5 ),
                                           Slope=signif(efit1$coef[[2]], 5),
                                           Estimate=signif(summary(efit1)$coef[2,1], 5),
                                           Std.Error=signif(summary(efit1)$coef[2,2], 5),
                                           t.value=signif(summary(efit1)$coef[2,3], 5),
                                           P=signif(summary(efit1)$coef[2,4], 5),
                                           N=dim(edata)[1]))
}

eN=erandom_tab$N

for (i in 1:10000){
  erows <- sample(1:dim(tab)[1], sample(eN,1), replace = F)
  etable <- data.frame(e_div=tab$e_div,fst=tab$fst)[erows,]#e_div values remain fixed and we shuffel the fst values to that:
  etable$fst=sample(etable$fst,dim(etable)[1],replace = T)
  efit <- lm(fst~e_div,data=etable)
  erandom_tab=rbind(erandom_tab,data.frame(group="random",
                                           R2=signif(summary(efit)$r.squared, 5),
                                           Adj_R2=signif(summary(efit)$adj.r.squared, 5),
                                           Intercept=signif(efit$coef[[1]],5 ),
                                           Slope=signif(efit$coef[[2]], 5),
                                           Estimate=signif(summary(efit)$coef[2,1], 5),
                                           Std.Error=signif(summary(efit)$coef[2,2], 5),
                                           t.value=signif(summary(efit)$coef[2,3], 5),
                                           P=signif(summary(efit)$coef[2,4], 5),
                                           N=dim(etable)[1]))
}




erandom_only=subset(erandom_tab,group=="random")
eoriginal=subset(erandom_tab,group!="random")


#P-value will be calculated as (# of random tvalues >= original tvalues) / (# of random tvalues < original tvalues)
ggplot(erandom_only,aes(x=t.value))+
  geom_histogram(binwidth=0.2)+
  geom_vline(xintercept=eoriginal$t.value)+
  ggtitle("e_div")

#please note that these plots and results might slightly change, depending on the random distribution of tvalues that was generated
adj.p.val = c()
for (i in 1:9){
  p =nrow(erandom_only[abs(erandom_only$t.value)>=eoriginal[i,8],])/nrow(erandom_only[abs(erandom_only$t.value)<eoriginal[i,8],])
  adj.p.val = c(adj.p.val,p)
}

#attach adjusted p values to table:
eoriginal$adj.p.val = adj.p.val





#validate randomization:
hist(erandom_only$P)
abline(v=eoriginal$P,col="red")
#

ggplot(erandom_only,aes(x=R2))+
  geom_histogram(binwidth=0.001)+
  geom_vline(xintercept=eoriginal$R2)

ggplot(erandom_only,aes(x=Adj_R2))+
  geom_histogram(binwidth=0.0005)+
  geom_vline(data=eoriginal, aes(xintercept=Adj_R2, col=group))+
  xlim(c(-0.003,0.065))+
  ggtitle("e_div")

ggplot(erandom_only,aes(x=Slope))+
  geom_histogram(binwidth=0.1)+
  geom_vline(data=eoriginal, aes(xintercept=Slope, col=group))+
  xlim(c(-1.5,8.5))+
  ggtitle("e_div")





write.csv2(eoriginal,"e_div_fst_regression_original_data.csv",row.names = T)
write.csv2(erandom_tab,"e_div_fst_regression_randomized_and_original.csv",row.names = T)





############### absolute variation explained by geographic distance (log_geo)
#similar to the approach above for evolutionary distance
#we are interested in the question: How much variation in fst is explained by log_geo?
#Thus regression is a good model to use
#data is not normally distributed, thats why we cannot trust the p value calculated for the regression
#with a bootstrap randomization test we create our own distribution
#randomization tests (10,000x bootstrap)

#create the groups you want to look at:
tab$type_introvect <- paste(tab$type,tab$introduction_vector,sep="_")

#there are na values in log_geo
tab2 <- na.omit(tab)


gfit_original <- lm(fst~log_geo,data=tab2)

#create table with original and randomized values:
grandom_tab=data.frame(group="All",
                       R2=signif(summary(gfit_original)$r.squared, 5),
                       Adj_R2=signif(summary(gfit_original)$adj.r.squared, 5),
                       Intercept=signif(gfit_original$coef[[1]],5 ),
                       Slope=signif(gfit_original$coef[[2]], 5),
                       Estimate=signif(summary(gfit_original)$coef[2,1], 5),
                       Std.Error=signif(summary(gfit_original)$coef[2,2], 5),
                       t.value=signif(summary(gfit_original)$coef[2,3], 5),
                       P=signif(summary(gfit_original)$coef[2,4], 5),
                       N=dim(tab2)[1])

###add the values for the 4 groups (introduction*type)
for (i in unique(tab2$type_introvect)){
  gdata = subset(tab2[tab2$type_introvect==i,])
  gfit1 <- lm(fst~log_geo,data=gdata)
  grandom_tab=rbind(grandom_tab,data.frame(group=i,
                                           R2=signif(summary(gfit1)$r.squared, 5),
                                           Adj_R2=signif(summary(gfit1)$adj.r.squared, 5),
                                           Intercept=signif(gfit1$coef[[1]],5 ),
                                           Slope=signif(gfit1$coef[[2]], 5),
                                           Estimate=signif(summary(gfit1)$coef[2,1], 5),
                                           Std.Error=signif(summary(gfit1)$coef[2,2], 5),
                                           t.value=signif(summary(gfit1)$coef[2,3], 5),
                                           P=signif(summary(gfit1)$coef[2,4], 5),
                                           N=dim(gdata)[1]))
}




gN=grandom_tab$N

for (i in 1:10000){
  grows <- sample(1:dim(tab2)[1], sample(gN,1), replace = F)
  gtable <- data.frame(log_geo=tab2$log_geo,fst=tab2$fst)[grows,]#e_div values remain fixed and we shuffel the fst values to that:
  gtable$fst=sample(gtable$fst,dim(gtable)[1],replace = T)
  gfit <- lm(fst~log_geo,data=gtable)
  grandom_tab=rbind(grandom_tab,data.frame(group="random",
                                           R2=signif(summary(gfit)$r.squared, 5),
                                           Adj_R2=signif(summary(gfit)$adj.r.squared, 5),
                                           Intercept=signif(gfit$coef[[1]],5 ),
                                           Slope=signif(gfit$coef[[2]], 5),
                                           Estimate=signif(summary(gfit)$coef[2,1], 5),
                                           Std.Error=signif(summary(gfit)$coef[2,2], 5),
                                           t.value=signif(summary(gfit)$coef[2,3], 5),
                                           P=signif(summary(gfit)$coef[2,4], 5),
                                           N=dim(gtable)[1]))
}



grandom_only=subset(grandom_tab,group=="random")
goriginal=subset(grandom_tab,group!="random")



#P-value will be calculated as (# of random tvalues >= original tvalues) / (# of random tvalues < original tvalues)
ggplot(grandom_only,aes(x=t.value))+
  geom_histogram(binwidth=0.2)+
  geom_vline(xintercept=goriginal$t.value)+
  ggtitle("log_geo")

#please note that these plots and results might slightly change, depending on the random distribution of tvalues that was generated
adj.p.val = c()
for (i in 1:9){
  p =nrow(grandom_only[abs(grandom_only$t.value)>=goriginal[i,8],])/nrow(grandom_only[abs(grandom_only$t.value)<goriginal[i,8],])
  adj.p.val = c(adj.p.val,p)
}

#attach adjusted p values to table:
goriginal$adj.p.val = adj.p.val





#validate randomization:
hist(grandom_only$P)
abline(v=goriginal$P,col="red")
#

ggplot(grandom_only,aes(x=R2))+
  geom_histogram(binwidth=0.001)+
  geom_vline(xintercept=goriginal$R2)

ggplot(grandom_only,aes(x=Adj_R2))+
  geom_histogram(binwidth=0.0005)+
  geom_vline(data=goriginal, aes(xintercept=Adj_R2, col=group))+
  xlim(c(-0.003,0.065))+
  ggtitle("log_geo")

ggplot(grandom_only,aes(x=Slope))+
  geom_histogram(binwidth=0.1)+
  geom_vline(data=goriginal, aes(xintercept=Slope, col=group))+
  xlim(c(-1.5,8.5))+
  ggtitle("log_geo")







write.csv2(goriginal,"log_geo_fst_regression_original_data.csv",row.names = T)
write.csv2(grandom_tab,"log_geo_fst_regression_randomized_and_original.csv",row.names = T)








####### relative proportion of variation in FST explained by evolutionary and geographic distance (IC)
#overall dataset:
IC_mod <- hier.part(tab2$fst,tab2[,c(5,6)], fam = "quasibinomial", gof =  "RMSPE")
#
#        ind.exp.var
#e_div      88.87395
#log_geo    11.12605
#-----> Table 4

#eight subgroups of interest:
tab2$type_introvect <- paste(tab2$type,tab2$introduction_vector,sep="_")

#to combine all value is one table, setup a table with the IC values for th eoverall dataset:
IC_tab=data.frame(group="All",
                  IC_ediv=IC_mod$I.perc[1,],
                  IC_geo=IC_mod$I.perc[2,],
                  N=dim(tab2)[1])

for (i in unique(tab2$type_introvect)){
  data = subset(tab2[tab2$type_introvect==i,])
  IC_mod <- hier.part(data$fst,data[,c(5,6)], fam = "quasibinomial", gof =  "RMSPE")
  IC_tab=rbind(IC_tab,data.frame(group=i,
                                 IC_ediv=IC_mod$I.perc[1,],
                                 IC_geo=IC_mod$I.perc[2,],
                                 N=dim(data)[1]))
  
}


#Next, test whether the subgroup IC values deviate for the overall group values using a chisquare test:
#in this case all will be our expected value and the subgroup values our observed values
#check for overall differences accross all values:
chisq.p <-chisq.test(IC_tab[,2:3])$p.value

for (i in 2:5){
  chisq.p <- c(chisq.p, chisq.test(IC_tab[c(1,i),2:3])$p.value)
}

#correct pvalues for multiple testing:
IC_tab$chisq.p.fdr=round(p.adjust(chisq.p, method = "fdr", n = length(chisq.p)),3)
#-----> the resulting IC_tab is part of Table 4




#write.csv2(IC_tab,"IC_values_all_groups.csv",row.names = F)
