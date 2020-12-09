#load packages
library("ggplot2")
library("cowplot")
library("lme4")
library("car")
library("multcomp")

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
ggplot(tab,aes(x=species,y=fst))+
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
#-----> Table 1 in the paper


#model validation:
plot(ModLog)
plot(residuals(ModLog)~tab2$lifespan_zscores)
#looks okay


ggplot(tab2,aes(x=lifespan_zscores,y=log10(m_geo)))+
  geom_point()













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
#-----> Table 2 in the paper





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





####we find a single_multiple_introduction main effect over the entire dataset (native and non-native)
# However, this factor should only affect the non-native populations and not the native populations
# Here we test if the single_multiple_introduction main effect is present in the native / non-native range sub-dataset
non.nat <- subset(tab,type=="non-native")
nat <- subset(tab,type=="native")

Mod2nn = glmer(fst ~ single_multiple_introductions+
                 (1|gene/species), family=binomial, data=non.nat, 
               control=glmerControl(optimizer="nloptwrap", optCtrl=list(xtol_abs=1e-12,ftol_abs=1e-12)))

summary(Mod2nn)
Anova(Mod2nn,type=2)
#                               Chisq Df Pr(>Chisq)  
#single_multiple_introductions 5.3175  1    0.02111 *

Mod2n = glmer(fst ~ single_multiple_introductions+
                (1|gene/species), family=binomial, data=nat, 
              control=glmerControl(optimizer="nloptwrap", optCtrl=list(xtol_abs=1e-12,ftol_abs=1e-12)))

summary(Mod2n)
Anova(Mod2n,type=2)
#                               Chisq Df Pr(>Chisq)
#single_multiple_introductions 1.8115  1     0.1783

#----> we can conclude:
# Comparing the number of introduction events, pairwise FST was significantly lower in populations 
#originating from a single introduction event compared to those from multiple introductions (Chisq1=4.97, P=0.026). 
#Further tests within the non-native and native range show that this effect holds true for the populations 
#living in their non-native range (Chisq1=5.318, P=0.021), but not for populations living 
#in their native range (Chisq1=1.812, P=0.178). (see publication for further details)















