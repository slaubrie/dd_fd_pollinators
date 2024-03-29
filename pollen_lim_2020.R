# last updated 21 March 2021 
# analysis to look at if bagging plants affects their success rate in four focal plants in perenjori 
# still figuring out github..... :| 

# libraries
require(lme4)
require(emmeans)
require(ggplot2)
require(glmmTMB)

# data 
dat<-read.csv('Bagged_Plants2020.csv', header=T)

# complete cases
df<-dat[which(complete.cases(dat)==T),]

# add failures in 
df$fails<-(df$attempts)-(df$successes)

# add in id 
df$id<-seq(1:nrow(df))

# two column integer matrix with successes then fails 
# for more info do: ?binomial
try<-cbind(df$successes, df$fails)

# model
mod<-glmmTMB(cbind(successes,fails)~species*treatment+(1|id), family='binomial', data=df)

# get data 
estimates<-emmip(mod, species~treatment, CI=T, type='response', plotit=F)

# plot 
plot(10,10)
labels<-c(arca="ARCA", vero="GORO", tror="TROR", trcy="TRCY")
ggplot(estimates, aes(x=xvar, y=yvar, group=species,col=species, shape=species))+
  geom_line(lwd=1.2)+
  geom_pointrange(aes(ymin=LCL, ymax=UCL), size=1)+
  theme_bw()+theme(text=element_text(size=20), 
                   axis.text.x=element_text(size=15), 
                   axis.text.y=element_text(size=15))+
  labs(col='Species', shape='Species')+
  scale_color_manual(labels=labels, values=c('#002E09', "#E7B800",  "#00AFBB", "#FC4E07"))+
  scale_shape_manual(labels=labels, values=c(16, 15, 3, 17))+
  xlab('Treatment')+
  ylab('Seed Success')

# get pairs comparison on the species two different bagging treatments
mod.emm<-emmeans(mod, ~treatment|species)
tab<-pairs(mod.emm)

# can also compare species success 
mod.emm2<-emmeans(mod, ~species|treatment)
tab2<-pairs(mod.emm2)

# prob of setting seeds
emmeans(mod, ~species|treatment, type="response")
