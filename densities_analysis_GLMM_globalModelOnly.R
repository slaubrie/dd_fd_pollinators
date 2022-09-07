# looking at density and frequency dependence in datasets from experiment run in 2020
# using one global GLMM model
# created 23 September 2021
# last updated 6 April 2022

# packages 
require(ggplot2)
require(gamm4)
require(dplyr)
require(visreg)
require(ggpubr)
require(mgcViz)
require(MuMIn) # this is required to get AICc
require(performance) #for getting r-squared values: https://easystats.github.io/performance/reference/r2_nakagawa.html
require(lme4)
require(MASS)
## this is how i installed glmmTMB remotes::install_github("glmmTMB/glmmTMB/glmmTMB")
require(glmmTMB)
require(stringr)
require(jtools)
require(sjPlot)

# get data 

# set working directory (for aub)
setwd('/Users/aubrie/Dropbox/UQ/WesternAustraliaProjects/PFD2020/DataCode/')

# data 
dat<-read.csv('pfdata_2020.csv', header=T)
data<-dat[,c("Species","Trim.Treatment2","rep","Fate","viable_sd","nviable_sd","ring_consp","ring_hetero","total_consp","total_hetero")]

#rep as factor
data$rep<-as.factor(data$rep) # random effect poss

# make the class of viable seed numeric, wtf??????????
data$viable_sd<-as.numeric(as.character(data$viable_sd))
data$nviable_sd<-as.numeric(as.character(data$nviable_sd))

# full rows
data<-data[which(complete.cases(data)==T),]

# rename trim treatment
names(data)[names(data)=='Trim.Treatment2']<-'Trim.Treatment'
data$Trim.Treatment<-as.factor(data$Trim.Treatment)

# log transform the densities
data$log_ringcon<-log1p(data$ring_consp)
data$log_ringhet<-log1p(data$ring_hetero)

# can look at frequencies too
data$freq.range<-data$ring_hetero/(rowSums(data[,7:8]))

# mean and standard deviation of frequencies 
mean(na.omit(data$freq.range))

sd(na.omit(data$freq.range))

### fitting some gams with gamm4 and then plotting with ggplot

## split by species
spl<-split(data, data$Species, drop=T)
trcy<-spl$TRCY
tror<-spl$TROR
vero<-spl$VERO
arca<-spl$ARCA

### some notes on glm analysis

# the models using poisson are overdispersed according to every check i use, so i'm using glmmTMB with nbinom2

## stealing code from ben bolker to determine if i need neg binom vs poisson http://bbolker.github.io/mixedmodels-misc/glmmFAQ.html#overdispersion
overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}

# could also do a LRT to compare the fit of the poisson and the negative binomial, because poisson is nested in negative binomial :-)

################## ################## #########
################## TRCY data ################## 
################## ################## #########

### model: 15 cm ring

## model comparison
mod.trcy.ring<-glmmTMB(viable_sd~poly(log_ringcon,2)*Trim.Treatment+poly(log_ringhet,2)*Trim.Treatment+(1|rep), data=trcy, family=nbinom2)

# summary
summary(mod.trcy.ring)
r2_nakagawa(mod.trcy.ring)

## predict

# generate data to predict model over 
trcy.newd.ring<-expand.grid(log_ringcon=seq(min(trcy$log_ringcon), max(trcy$log_ringcon), length.out=nrow(trcy)), log_ringhet=seq(min(trcy$log_ringhet), max(trcy$log_ringhet), length.out=nrow(trcy)), Trim.Treatment=c("CUT","NOCUT"), rep=as.factor(c(7:12))) 

# predict values of seed set 
trcy.pred.ring<-predict(mod.trcy.ring, newdata=trcy.newd.ring, type='response', re.form=NA)

# dataframe
trcy.newd.ring<-cbind(trcy.newd.ring, trcy.pred.ring)

################## ################## #########
################## TROR data ################## 
################## ################## #########

# remove insane outliers (there are three of them...)
tror<-tror[which(tror$viable_sd<170),]

### model: 15cm ring

mod.tror.ring<-glmmTMB(viable_sd~poly(log_ringcon,2)*Trim.Treatment+poly(log_ringhet,2)*Trim.Treatment+(1|rep), data=tror, family=nbinom2)

r2_nakagawa(mod.tror.ring)

# summary
summary(mod.tror.ring)

# generate data to predict model over 
tror.newd.ring<-expand.grid(log_ringcon=seq(min(tror$log_ringcon), max(tror$log_ringcon), length.out=nrow(tror)), log_ringhet=seq(min(tror$log_ringhet), max(tror$log_ringhet), length.out=nrow(tror)), Trim.Treatment=c("CUT","NOCUT"), rep=as.factor(c(13:18))) 

# predict values of seed set
tror.pred.ring<-predict(mod.tror.ring, newdata=tror.newd.ring, type='response', re.form=NA)

# dataframe
tror.newd.ring<-cbind(tror.newd.ring, tror.pred.ring)

################## ################## #########
################## VERO data ################## 
################## ################## #########

### model  15cm

mod.vero.ring<-glmmTMB(viable_sd~poly(log_ringcon,2)*Trim.Treatment+poly(log_ringhet,2)*Trim.Treatment+(1|rep), data=vero, family=nbinom2)

# summary
summary(mod.vero.ring)
r2_nakagawa(mod.vero.ring)

# generate data to predict model over 
vero.newd.ring<-expand.grid(log_ringcon=seq(min(vero$log_ringcon), max(vero$log_ringcon), length.out=nrow(vero)), log_ringhet=seq(min(vero$log_ringhet), max(vero$log_ringhet), length.out=nrow(vero)), Trim.Treatment=c("CUT","NOCUT"), rep=as.factor(c(19:24))) 

# predict values of seed set
vero.pred.ring<-predict(mod.vero.ring, newdata=vero.newd.ring, type='response', re.form=NA)

# dataframe
vero.newd.ring<-cbind(vero.newd.ring,vero.pred.ring)

################## ################## #########
################## ARCA data ################## 
################## ################## #########

# model 15 cm 

mod.arca.ring<-glmmTMB(viable_sd~poly(log_ringcon,2)*Trim.Treatment+poly(log_ringhet,2)*Trim.Treatment+(1|rep), data=arca, family=nbinom2)

# summary
summary(mod.arca.ring) 
r2_nakagawa(mod.arca.ring)

# generate data to predict model over 
arca.newd.ring<-expand.grid(log_ringcon=seq(min(arca$log_ringcon), max(arca$log_ringcon), length.out=nrow(arca)), log_ringhet=seq(min(arca$log_ringhet), max(arca$log_ringhet), length.out=nrow(arca)), Trim.Treatment=c("CUT","NOCUT"), rep=as.factor(c(1:6))) 

# predict values of seed set
arca.pred.ring<-predict(mod.arca.ring, newdata=arca.newd.ring, type='response', re.form=NA)

# dataframe
arca.newd.ring<-cbind(arca.newd.ring,arca.pred.ring)

###############################
######## summaries ############
###############################

## TRCY
summary(mod.trcy.ring)
r2_nakagawa(mod.trcy.ring)

## TROR
summary(mod.tror.ring)
r2_nakagawa(mod.tror.ring)

## VERO
summary(mod.vero.ring)
r2_nakagawa(mod.vero.ring)

## ARCA
summary(mod.arca.ring)
r2_nakagawa(mod.arca.ring)

###############################
############ PLOTS ############
###############################

## *~**~*~** CoLoUrS **~*~**~*
john.ramp.cols<- c("steelblue4", "lightsteelblue1", "yellow", "orange", "red")

##### ##### ##### ##### 
##### SMALL SCALE ##### 
##### ##### ##### ##### 

### TRCY 15cm 
b<-ggplot(trcy.newd.ring, 
          aes(log_ringcon, log_ringhet))+
  geom_tile(aes(fill=trcy.pred.ring))+
  stat_contour(bins=12,aes(log_ringcon, log_ringhet,z=trcy.pred.ring), 
               color="black", size=0.3, alpha=0.3)+
  scale_fill_gradientn(colours=john.ramp.cols, name="TRCY\nfecundity")+
  labs(x="", y="")+
  theme(text=element_text(size=20), 
        axis.text.x=element_text(size=15), 
        axis.text.y=element_text(size=15))+
  geom_jitter(aes(x=log_ringcon, y=log_ringhet, size=viable_sd),
              trcy, width=0.1, height=0.1, 
              alpha=0.5, 
              color='black', 
              pch=21,
              show.legend = F)
              
  #facet_wrap(vars(Trim.Treatment))
b

### TROR 15cm 
d<-ggplot(tror.newd.ring, 
          aes(log_ringcon, log_ringhet))+
  geom_tile(aes(fill= tror.pred.ring))+
  stat_contour(bins=12,aes(log_ringcon, 
                           log_ringhet,z=tror.pred.ring), 
               color="black", size=0.3, alpha=0.3)+
  scale_fill_gradientn(colours=john.ramp.cols, 
                       name='TROR\nfecundity')+
  labs(x="", y="")+theme(text=element_text(size=20), 
                         axis.text.x=element_text(size=15), 
                         axis.text.y=element_text(size=15))+
  geom_jitter(aes(x=log_ringcon, y=log_ringhet, size=viable_sd),
              tror, width=0.1, height=0.1, 
              alpha=0.5, color='black', pch=21,
              show.legend = F)

  #facet_wrap(vars(Trim.Treatment))
d

####### VERO 15cm
f<-ggplot(vero.newd.ring, 
          aes(log_ringcon, log_ringhet))+
  geom_tile(aes(fill=vero.pred.ring))+
  stat_contour(bins=12,aes(log_ringcon, log_ringhet,z=vero.pred.ring), 
               color="black", size=0.3, alpha=0.3)+
  scale_fill_gradientn(colours=john.ramp.cols, 
                       name="VERO\nfecundity")+
  labs(x="", y="")+theme(text=element_text(size=20), 
                         axis.text.x=element_text(size=15), 
                         axis.text.y=element_text(size=15))+
  facet_wrap(vars(Trim.Treatment))+
  geom_jitter(aes(x=log_ringcon, y=log_ringhet, size=viable_sd),
              vero, width=0.1, height=0.1,
              alpha=0.5, 
              color='black', 
              pch=21,
              show.legend = F)
f

####### ARCA 15cm
h<-ggplot(arca.newd.ring, 
          aes(log_ringcon, log_ringhet))+
  geom_tile(aes(fill=arca.pred.ring))+
  stat_contour(bins=12,aes(log_ringcon, log_ringhet,z=arca.pred.ring), 
               color="black", size=0.3, alpha=0.3)+
  scale_fill_gradientn(colours=john.ramp.cols, 
                       name="ARCA\nfecundity")+
  labs(x="", y="")+theme(text=element_text(size=20), 
                         axis.text.x=element_text(size=15), 
                         axis.text.y=element_text(size=15))+
  geom_jitter(aes(x=log_ringcon, y=log_ringhet, size=viable_sd),
              arca, width=0.1, height=0.1, 
              alpha=0.5, 
              color='black', 
              pch=21,
              show.legend = F)

  #facet_wrap(vars(Trim.Treatment))
h

ring<-ggarrange(h,b,d,f, nrow=2, ncol=2)
ring.test<-annotate_figure(ring, bottom=text_grob('log(Conspecific Density)+1', size=20), left=text_grob('log(Heterospecific Density)+1',size=20, rot=90))
ring.test
#ggsave(file='ring_bestfit.png', plot=ring.test, width=10, height=8, units="in")

#######################
##### PLOT COEFFS #####
#######################
coeff_plot_2020<-plot_models(mod.arca.ring, mod.trcy.ring, mod.tror.ring, mod.vero.ring, 
                            vline.color='gray', 
                             colors=c( '#002E09', "#FC4E07", "#E7B800","#00AFBB"), 
                             m.labels = c("ARCA", "TRCY", "TROR", "VERO"), axis.title=("Estimate"), legend.title = ("Model"), transform=NULL)+
  theme_bw()+
  theme(axis.text.x=element_text(size=17), 
        axis.text.y=element_text(size=15), 
        axis.title.y=element_text(size=15), 
        axis.title.x=element_text(size=20), 
        legend.title=element_text(size=20), 
        legend.text=element_text(size=15))
#coeff_plot_2020
#ggsave(file='coeff_plot_2020.png', plot=coeff_plot_2020, width=8, height=6, units="in")

#######################
##### COEFF TABLE #####
#######################

sjPlot::tab_model(mod.arca.ring, mod.trcy.ring, mod.tror.ring, mod.vero.ring, 
                  transform=NULL, 
                  dv.labels=c("ARCA","TRCY","TROR","VERO"),string.pred = "Coeffcient",
                  string.ci = "Conf. Int (95%)",
                  string.p = "P-Value")
