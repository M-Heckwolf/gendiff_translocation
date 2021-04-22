#Script for "Fresh start after rough rides: understanding patterns of genetic differentiation upon human-mediated translocations"
#Melanie Heckwolf*, Teófilo Morim*, Francesca Riccioli, Miguel Baltazar-Soares (*shared first authors)

#load packages
library("ggplot2")
library("cowplot")
library("lme4")
library("car")
library("multcomp")
library("dplyr")

#set theme for plots
theme_set(theme_cowplot())

#load table
tab <- read.csv2("\\\\path_to_file/global_data.csv")

#explore table
str(tab)
summary(tab)





#
#
#
#
############ STATISTICAL MODELS #################################################
hist(tab$fst) #distribution binomial
#run binomial model with fst as response variable:

#random factors: gene / species
#different papers used differend species / genes and they differ in their fst values
#we nest species in gene, since each species was only tested for one gene and we have fewer genes then species
ggplot(tab,aes(x=species,y=fst,fill=parental_care))+
  geom_boxplot()

ggplot(tab,aes(x=species,y=fst,col=single_multiple_introductions))+
  geom_boxplot()

ggplot(tab,aes(x=gene,y=fst))+
  geom_boxplot()

#fixed factors: 
# type, parental care, intro vector, SM intro, feeding,
# dispersal, lifespan
# interactions: type x dispersal, type x intro vector, type x SM intro

#we did not include e_div, log_geo since we will focus on these variables later
#we did not include group since it is a very artificial parameter with sometimes no 
#    species replication within one group so it would not be enough data to
#    extrapolate these findings to a group level
#we did not include time since introduction as this variable is only meaningful for the non-native group, but
#    we did also no include the interaction between type and time since introduction
#    since the data for each type along the time axis is very unevenly distributed






#------------------> Model 1: differences in sampling area between the groups:

#to ensure that differences in fst are not caused by differences in sampling are between the groups
#we run the full and the reduced model with the same fixed and random factors but log_geo as a response variable:
#NA values in log_geo:
tab2=na.omit(tab)
tab2$km_geo=10**tab2$log_geo
tab2$m_geo=round(tab2$km_geo*1000)

#distribution of log_geo values is skewed (also after transformation)
#we run a model with gamma error distribution on log transformed distances in meters

ModLog = glmer(log10(m_geo) ~ type+parental_care+introduction_vector+single_multiple_introductions+
                 feeding+dispersal+lifespan_zscores+(1|gene/species), family=Gamma, data=tab2)

derivs1 <- ModLog@optinfo$derivs
sc_grad1 <- with(derivs1,solve(Hessian,gradient))
max(abs(sc_grad1)) #0.0003392479 should be <0.001, looks good
max(pmin(abs(sc_grad1),abs(derivs1$gradient)))
#Warning is not an issue, I tried also other optimzer and xtol/ftol values and they give the same result


summary(ModLog)
Anova(ModLog,type=2)
#                                Chisq Df Pr(>Chisq)    
#type                           1.4963  1  0.2212481    
#parental_care                  0.5213  1  0.4703023    
#introduction_vector            0.0018  1  0.9663668    
#single_multiple_introductions  0.4671  1  0.4943153    
#feeding                        6.0669  3  0.1083986    
#dispersal                      0.0233  1  0.8785712    
#lifespan_zscores              12.5790  1  0.0003901 ***  
#-----> Table S5 in the Manuscript (supplementary table)


#model validation:
plot(ModLog)
plot(residuals(ModLog)~tab2$lifespan_zscores)
#looks okay


ggplot(tab2,aes(x=lifespan_zscores,y=log10(m_geo)))+
  geom_point()





################### we will further validate our results by:
################### (1.1) running the same model on 1000 rarefied datasets, to account for uneven sample size
################### (1.2) running an ANOVA on median geographic distance per species and range (native, non-native)
################### (1.3) repeating the GLMM analysis excluding samples that were present in only one of the two ranges (native, non-native)

####### (1.1)
# the 20 studies / species are represented at a different frequency within this dataset (from 3 to 990 fst per species and range)
# to ensure the results are not biased through uneven sample sizes between studies we ran the analysis again with rarefied datasets (6 per species and range, representing the smallest 5% of the groups)
tab2$species_type <- paste(tab2$species, tab2$type, sep="_")
table(table(tab2$species_type))# 3 to 990 comparisons per species and range


#save original read statistics in a table and add rarefied dataset analyses:
Mod.Geo <- data.frame()
for (i in 1:3){
  row <- cbind(data.frame(data="original",
                              run=NA,
                              variable=rownames(t(Anova(ModLog,type=2)[i]))),
                   t(Anova(ModLog,type=2)[i]))
  Mod.Geo <- rbind(Mod.Geo,row)
}

# run rarefaction analysis:
a=0
for (i in 1:1000){
  a=a+1 #to add run number later and keep track of run
  print(a)
  rarefied <- tab2[tab2$species_type!="Microcosmus squamiger_native",] %>% group_by(species_type) %>% sample_n(6) #subsample dataset to 3 comparisons per group
  rarefied <- rbind(rarefied,tab2[tab2$species_type=="Microcosmus squamiger_native",])
  ModLog_rare = glmer(log10(m_geo) ~ type+parental_care+introduction_vector+single_multiple_introductions+
                        feeding+dispersal+lifespan_zscores+(1|gene/species), family=Gamma, data=rarefied)
  
  for (i in c(1,3)){
    row <- cbind(data.frame(data="random",
                            run=a,
                            variable=rownames(t(Anova(ModLog_rare,type=2)[i]))),
                 t(Anova(ModLog_rare,type=2)[i]))
    Mod.Geo <- rbind(Mod.Geo,row)
  }
  
}

#write.csv2(Mod.Geo,"Rarefied_LogGeo_GLMMgamma_7fcts_N6.csv",row.names=FALSE)
# --> Table S6 in the Manuscript (supplementary table)






####### (1.2)
# run an ANOVA on median geogrpahic distances to account for the fact that in pariwise geographic distance measures, single populations are used multiple times.
# This fact can create pseudo-replication issues, thus here we verify that our results are robust against that.
tab2 <- tab2 %>% group_by(species_type) %>% mutate(median_geo=median(m_geo))
tab2.means <- tab2[,c("type","parental_care","introduction_vector",
                      "single_multiple_introductions","feeding","median_geo",
                      "dispersal","lifespan_zscores","species_type")] %>% distinct()


ModLog2II = lm(log10(median_geo) ~ type+parental_care+introduction_vector+single_multiple_introductions+
                 feeding+dispersal+lifespan_zscores, tab2.means)

summary(ModLog2II)
anova(ModLog2II)
#                              Df  Sum Sq Mean Sq F value   Pr(>F)   
#type                           1 0.01030 0.01030  0.1159 0.736928   
#parental_care                  1 0.19017 0.19017  2.1401 0.158296   
#introduction_vector            1 0.15503 0.15503  1.7446 0.200765   
#single_multiple_introductions  1 0.36078 0.36078  4.0601 0.056888 . 
#feeding                        3 0.47046 0.15682  1.7648 0.184699   
#dispersal                      1 0.38400 0.38400  4.3214 0.050081 . 
#lifespan_zscores               1 1.27558 1.27558 14.3551 0.001075 **
#Residuals                     21 1.86603 0.08886                   

#--> Table 1 in the Manuscript








####### (1.3)
#repeat the original GLMM only with the 11 species we have inforation about in their native AND non-native range:
tab2.pairs <- subset(tab2, species %in% c("Barbus barbus","Coleophora deauratella","Cyprinella lutrensis",
                                          "Gambusia holbrooki","Hemimysis anomala","Metrioptera roeselii",
                                          "Microcosmus squamiger","Mnemiopsis leidyi","Pacifastacus leniusculus",
                                          "Pseudorasbora parva","Podarcis siculus"))

#we run a model with gamma error distribution on log transformed distances in meters
ModLog.pairs = glmer(log10(m_geo) ~ type+parental_care+introduction_vector+single_multiple_introductions+
                 feeding+dispersal+lifespan_zscores+(1|gene/species), family=Gamma, data=tab2.pairs)

summary(ModLog.pairs)
Anova(ModLog.pairs,type=2)
#                                Chisq Df Pr(>Chisq)    
#type                           2.8316  1  0.0924261 .  
#parental_care                  0.7670  1  0.3811520    
#introduction_vector            0.0954  1  0.7574342    
#single_multiple_introductions  0.4602  1  0.4975389    
#feeding                        0.1698  1  0.6802536    
#dispersal                      2.3681  1  0.1238349    
#lifespan_zscores              13.6759  1  0.0002172 ***

# --> Table S7 in the Manuscript (supplementary table)












#------------------> Model 2: genetic differentiation between the groups:
#note that we do not use the counts for 0/1 in a cbind command but the fst values (from 0-1)
#    when coding it this way this warning message will appear and is normal:
#    1: In eval(family$initialize, rho) :
#    non-integer #successes in a binomial glm!

Mod2 = glmer(fst ~ type+parental_care+introduction_vector+single_multiple_introductions+
               feeding+dispersal+
               type:dispersal+type:introduction_vector+type:single_multiple_introductions+
               (1|gene/species), family=binomial, data=tab, 
             control=glmerControl(optimizer="nloptwrap", optCtrl=list(xtol_abs=1e-12,ftol_abs=1e-12)))


#
#Warning:
#    boundary (singular) fit: see ?isSingular
#test if it is singular with a tolerance of 1e-5
#-05 a widely accepted tolerance and the default for some functions (eg isSingular) 
#while -04 is the default tolerance within glmer
isSingular(Mod2, tol = 1e-5)#FALSE with this tolerance
#model is converging at the boundary but within our tolerance range
#this means we can ignore the warning

summary(Mod2)
Anova(Mod2,type=3)#here the interaction is taken into account, thus the values for the single factors with non sign. interactions differ (model selection should be done here)
#                                      Chisq Df Pr(>Chisq)    
#(Intercept)                          0.0223  1    0.88134    
#type                               121.8048  1  < 2.2e-16 ***
#parental_care                        0.9149  1    0.33883    
#introduction_vector                  0.3317  1    0.56468    
#single_multiple_introductions        0.0118  1    0.91343    
#feeding                              2.1625  3    0.53938    
#dispersal                            2.6267  1    0.10508    
#type:dispersal                       4.2300  1    0.03972 *  
#type:introduction_vector            41.8632  1  9.789e-11 ***
#type:single_multiple_introductions   0.0012  1    0.97289 
#-----> same as Table 2 but before model selection


#model validation:
plot(Mod2)
boxplot(residuals(Mod2)~tab$type)
boxplot(residuals(Mod2)~tab$dispersal)
boxplot(residuals(Mod2)~tab$introduction_vector)


########### model selection: remove type:single_multiple_introductions

Mod2b = glmer(fst ~ type+parental_care+introduction_vector+single_multiple_introductions+
                feeding+dispersal+
                type:dispersal+type:introduction_vector+
                (1|gene/species), family=binomial, data=tab, 
              control=glmerControl(optimizer="nloptwrap", optCtrl=list(xtol_abs=1e-12,ftol_abs=1e-12)))


#
#Warning:
#    boundary (singular) fit: see ?isSingular
#test if it is singular with a tolerance of 1e-5
#-05 a widely accepted tolerance and the default for some functions (eg isSingular) 
#while -04 is the default tolerance within glmer
isSingular(Mod2b, tol = 1e-7)#FALSE with this tolerance
#model is converging at the boundary but within our tolerance range
#this means we can ignore the warning


#model validation:
plot(Mod2b)
boxplot(residuals(Mod2b)~tab$type)
boxplot(residuals(Mod2b)~tab$dispersal)
boxplot(residuals(Mod2b)~tab$introduction_vector)

#compare models:
anova(Mod2,Mod2b)#Mod2b is the better fit model

#results
summary(Mod2b)
Anova(Mod2b,type=3)#here the interaction is taken into account, thus the values for the single factors with non sign. interactions differ (model selection should be done here)
#                                 Chisq Df Pr(>Chisq)    
#(Intercept)                     0.0212  1    0.88424    
#type                          121.7787  1  < 2.2e-16 ***
#parental_care                   0.9316  1    0.33444    
#introduction_vector             0.3319  1    0.56456    
#single_multiple_introductions   4.9731  1    0.02574 *  
#feeding                         2.2427  3    0.52359    
#dispersal                       2.6429  1    0.10401    
#type:dispersal                  4.2460  1    0.03934 *  
#type:introduction_vector       41.9848  1  9.199e-11 ***

#-----> Table S8 in the Manuscript (Supplementary table)






################### we will further validate our results by:
################### (1.1) running the same model on 1000 rarefied datasets, to account for uneven sample size
################### (1.2) running an ANOVA on median geographic distance per species and range (native, non-native)
################### (1.3) repeating the GLMM analysis excluding samples that were present in only one of the two ranges (native, non-native)

####### (1.1)
# the 20 studies / species are represented at a different frequency within this dataset (from 3 to 990 fst per species and range)
# to ensure the results are not biased through uneven sample sizes between studies we ran the analysis again with rarefied datasets (6 per species and range, representing the smallest 5% of the groups)
tab$species_type <- paste(tab$species, tab$type, sep="_")
table(tab$species_type)# 3 to 990 comparisons per species and range


#save original read statistics in a table and add rarefied dataset analyses:
Mod.gendiff <- data.frame()
for (i in 1:3){
  row <- cbind(data.frame(data="original",
                          run=NA,
                          variable=rownames(t(Anova(Mod2b,type=3)[i]))),
               t(Anova(Mod2b,type=3)[i]))
  Mod.gendiff <- rbind(Mod.gendiff,row)
}

# run rarefaction analysis:
a=0
for (i in 1:1000){
  a=a+1 #to add run number later and keep track of run
  print(a)
  rarefied <- tab[tab$species_type!="Microcosmus squamiger_native",] %>% group_by(species_type) %>% sample_n(6) #subsample dataset to 3 comparisons per group
  rarefied <- rbind(rarefied,tab[tab$species_type=="Microcosmus squamiger_native",])
  Mod2b_rare = glmer(fst ~ type+parental_care+introduction_vector+single_multiple_introductions+
                            feeding+dispersal+
                            type:dispersal+type:introduction_vector+
                            (1|gene/species), family=binomial, data=rarefied, 
                     control=glmerControl(optimizer="nloptwrap", optCtrl=list(xtol_abs=1e-12,ftol_abs=1e-12)))
  
  for (i in c(1,3)){
    row <- cbind(data.frame(data="random",
                            run=a,
                            variable=rownames(t(Anova(Mod2b_rare,type=3)[i]))),
                 t(Anova(Mod2b_rare,type=3)[i]))
    Mod.gendiff <- rbind(Mod.gendiff,row)
  }
  
}



#write.csv2(Mod.gendiff2,"Rarefied_FST_GLMMbinom_7fcts_N3.csv",row.names=FALSE)
# --> Table S9 in the Manuscript (supplementary table)








####### (1.2)
# run an ANOVA on median geogrpahic distances to account for the fact that in pariwise geographic distance measures, single populations are used multiple times.
# This fact can create pseudo-replication issues, thus here we verify that our results are robust against that.
tab <- tab %>% group_by(species_type) %>% mutate(mean_fst=median(fst))
tab.means <- tab[,c("type","parental_care","introduction_vector",
                    "single_multiple_introductions","feeding","mean_fst",
                    "dispersal","species_type")] %>% distinct()

Mod2b_2II <- aov(mean_fst ~ type+parental_care+introduction_vector+single_multiple_introductions+
                   feeding+dispersal+type:dispersal+type:introduction_vector, tab.means)
Anova(Mod2b_2II,type=3)
#                               Sum Sq Df F value    Pr(>F)    
#(Intercept)                   1.57429  1 22.3616 0.0001285 ***
#type                          0.43589  1  6.1915 0.0217712 *  
#parental_care                 0.05257  1  0.7468 0.3977469    
#introduction_vector           0.07941  1  1.1280 0.3008642    
#single_multiple_introductions 0.26760  1  3.8011 0.0653799 .  
#feeding                       0.01956  3  0.0926 0.9632615    
#dispersal                     0.22363  1  3.1764 0.0898961 .  
#type:dispersal                0.00957  1  0.1359 0.7162956    
#type:introduction_vector      0.46913  1  6.6637 0.0178289 *  
#Residuals                     1.40803 20                 


#model selection:
Mod2b_2II2 <- aov(mean_fst ~ type+parental_care+introduction_vector+single_multiple_introductions+
                    feeding+dispersal+type:introduction_vector, tab.means)
Anova(Mod2b_2II2,type=3)
anova(Mod2b_2II,Mod2b_2II2)
# 

Mod2b_2II3 <- aov(mean_fst ~ type+parental_care+introduction_vector+single_multiple_introductions+
                    dispersal+type:introduction_vector, tab.means)
Anova(Mod2b_2II3,type=3)
anova(Mod2b_2II,Mod2b_2II2,Mod2b_2II3)
#

Mod2b_2II4 <- aov(mean_fst ~ type+introduction_vector+single_multiple_introductions+
                    dispersal+type:introduction_vector, tab.means)
Anova(Mod2b_2II4,type=3)
anova(Mod2b_2II,Mod2b_2II2,Mod2b_2II3,Mod2b_2II4)
#

Mod2b_2II5 <- aov(mean_fst ~ type+single_multiple_introductions+
                    dispersal+type:introduction_vector, tab.means)
Anova(Mod2b_2II5,type=3)
anova(Mod2b_2II,Mod2b_2II2,Mod2b_2II3,Mod2b_2II4,Mod2b_2II5)
# same as model4, use model 4 as final model:

##### final model:
Anova(Mod2b_2II4,type=3)
#                               Sum Sq Df F value    Pr(>F)    
#(Intercept)                   2.06528  1 34.5497 3.927e-06 ***
#type                          0.53834  1  9.0058  0.006024 ** 
#introduction_vector           0.11285  1  1.8878  0.181642    
#single_multiple_introductions 0.38452  1  6.4326  0.017827 *  
#dispersal                     0.25936  1  4.3387  0.047639 *  
#type:introduction_vector      0.63253  1 10.5814  0.003263 ** 
#Residuals                     1.49443 25

# --> Table 2 in the Manuscript






####### (1.3)
#repeat the original GLMM only with the 11 species we have inforation about in their native AND non-native range:
tab.pairs <- subset(tab, species %in% c("Barbus barbus","Coleophora deauratella","Cyprinella lutrensis",
                                          "Gambusia holbrooki","Hemimysis anomala","Metrioptera roeselii",
                                          "Microcosmus squamiger","Mnemiopsis leidyi","Pacifastacus leniusculus",
                                          "Pseudorasbora parva","Podarcis siculus"))

#rerun the original model with these species only:
Mod2b.pairs = glmer(fst ~ type+parental_care+introduction_vector+single_multiple_introductions+
                feeding+dispersal+
                type:dispersal+type:introduction_vector+
                (1|gene/species), family=binomial, data=tab.pairs, 
              control=glmerControl(optimizer="nloptwrap", optCtrl=list(xtol_abs=1e-12,ftol_abs=1e-12)))

summary(Mod2b.pairs)
Anova(Mod2b.pairs,type=3)
#                                 Chisq Df Pr(>Chisq)    
#(Intercept)                     0.0959  1    0.75685    
#type                          119.2458  1  < 2.2e-16 ***
#parental_care                   0.0004  1    0.98322    
#introduction_vector             0.3826  1    0.53623    
#single_multiple_introductions   5.0837  1    0.02415 *  
#feeding                         0.0828  1    0.77350    
#dispersal                       3.3816  1    0.06593 .  
#type:dispersal                  2.8923  1    0.08900 .  
#type:introduction_vector       41.4824  1  1.189e-10 ***

# --> Table S10 in the Manuscript (supplementary table)








############## Post-hoc tests:

# recode interactions:
tab$Tintro <- interaction(tab$type,tab$introduction_vector,drop=T)
tab$Tdisp <- interaction(tab$type,tab$dispersal,drop=T)

#same model as before but different coding to enable post-hoc test:
Mod2ph = glmer(fst ~ Tdisp+type+parental_care+introduction_vector+single_multiple_introductions+
                 feeding+dispersal+
                 type:introduction_vector+
                 (1|gene/species), family=binomial, data=tab, 
               control=glmerControl(optimizer="nloptwrap", optCtrl=list(xtol_abs=1e-12,ftol_abs=1e-12)))

summary(glht(Mod2ph, mcp(Tdisp="Tukey")))
#                                            Estimate Std. Error z value Pr(>|z|)    
#non-native.active - native.active == 0       -1.6335     0.1481 -11.030   <0.001 ***
#native.passive - native.active == 0          -1.1274     0.7028  -1.604   0.3229    
#non-native.passive - native.active == 0      -2.0660     0.7384  -2.798   0.0199 *  
#native.passive - non-native.active == 0       0.5061     0.7162   0.707   0.8711    
#non-native.passive - non-native.active == 0  -0.4325     0.7288  -0.593   0.9185    
#non-native.passive - native.passive == 0     -0.9386     0.3542  -2.650   0.0302 *   
#-----> Table 3 in the paper

Mod2ph2 = glmer(fst ~ Tintro+parental_care+introduction_vector+single_multiple_introductions+
                  feeding+dispersal+type+
                  type:dispersal+
                  (1|gene/species), family=binomial, data=tab,
                control=glmerControl(optimizer="nloptwrap", optCtrl=list(xtol_abs=1e-12,ftol_abs=1e-12)))
#
summary(glht(Mod2ph2, mcp(Tintro="Tukey")))
#                                                       Estimate Std. Error z value Pr(>|z|)    
#non-native.intentional - native.intentional == 0       -1.63357    0.14802 -11.036   <1e-04 ***
#native.unintentional - native.intentional == 0         -0.34609    0.59926  -0.578   0.9275    
#non-native.unintentional - native.intentional == 0      0.01928    0.61836   0.031   1.0000    
#native.unintentional - non-native.intentional == 0      1.28747    0.60028   2.145   0.1129    
#non-native.unintentional - non-native.intentional == 0  1.65285    0.61379   2.693   0.0273 *  
#non-native.unintentional - native.unintentional == 0    0.36537    0.28310   1.291   0.5220    
#-----> Table 3 in the paper


