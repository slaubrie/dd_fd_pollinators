# Created in July 2021
# added to github October 2021
# analysis to look at if bagging plants affects their success rate in four focal plants in perenjori 

# libraries
require(lme4)
require(emmeans)
require(ggplot2)

# working directory
setwd('/Users/aubrie/Dropbox/UQ/WA.data/PFD2020/DataCode/')

# data 
dat<-read.csv('Bagged_Plants2020.csv', header=T)

# complete cases
df<-dat[which(complete.cases(dat)==T),]

# add failures in 
df$fails<-(df$attempts)-(df$successes)

# two column integer matrix with successes then fails 
# for more info do: ?binomial
try<-cbind(df$successes, df$fails)

# model
attach(df)
mod<-glm(try~species*treatment, family='binomial')

# get data 
estimates<-emmip(mod, species~treatment, CI=T, type='response', plotit=F)

# plot 
plot(10,10)
labels<-c(arca="ARCA", vero="VERO", tror="TROR", trcy="TRCY")
ggplot(estimates, aes(x=xvar, y=yvar, group=species,col=species))+geom_line()+geom_pointrange(aes(ymin=LCL, ymax=UCL))+theme_bw()+theme(text=element_text(size=20), axis.text.x=element_text(size=15), axis.text.y=element_text(size=15))+labs(col='Species')+scale_color_manual(labels=labels, values=c("#00AFBB", "#E7B800", "#FC4E07", '#002E09'))+xlab('Treatment')+ylab('Seed Success')

# get pairs comparison on the species two different bagging treatments
mod.emm<-emmeans(mod, ~treatment|species)
tab<-pairs(mod.emm)

# can also compare species success 
mod.emm2<-emmeans(mod, ~species|treatment)
tab2<-pairs(mod.emm2)

# prob of setting seeds
emmeans(mod, ~species|treatment, type="response")
