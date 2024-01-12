### doing vital rate analysis

# packages 
require(lme4)
require(emmeans)
require(pscl)
require(glmmTMB)
require(tidyr)
require(DHARMa)
require(ggplot2)
require(AICcmodavg)
require(ggpubr)

############## 2021 data #############

# read csv
dat<-read.csv('nplants_data_2021_git.csv', header=T)
dat1<-dat[which(dat$seeding_trt==1),]
dat1$physical_barrier<-as.factor(dat1$physical_barrier)
dat1$block<-as.factor(dat1$block)

############# treatment response: do analysis for count by species ############# 

dat2<-(dat1[c(1:4,7,11,15,22)]) # these are block, transect, initial, current_plot_type, ngoro, ntrcy, ntror, physical_barrier
head(dat2) 
countdat<-as.data.frame(dat2 %>% pivot_longer(c(ntror, ngoro, ntrcy)))
countdat$value<-as.numeric(ifelse(countdat$value>15, 15, countdat$value))
countmod<-glmmTMB(value~name*current_plot_type+(1|block), family="poisson", data=countdat)

# test for fit and zero inflation
 sim<-simulateResiduals(countmod)
 testZeroInflation(sim) # zero-inflated so need something else 
 plot(sim)

## going to do a hurdle model, which assumes a zero is only generated in one way 
# https://jsdajournal.springeropen.com/articles/10.1186/s40488-021-00121-4
# #https://stats.stackexchange.com/questions/81457/what-is-the-difference-between-zero-inflated-and-hurdle-models

# hurdle model 
# using examples as presented here: https://www.biorxiv.org/content/biorxiv/suppl/2017/05/01/132753.DC1/132753-2.pdf 
fit3<-glmmTMB(value~name*current_plot_type+(1|block), ziformula=~., family=nbinom2(), data=countdat)
 sim3<-simulateResiduals(fit3)
 plot(sim3)
 testDispersion(sim3) # looks nice ! 
 testZeroInflation(sim3) # looks nice ! 

summary(fit3)
emmip(fit3,~current_plot_type|name, type='response',CI=T)

# visualize
est<-emmeans(fit3,~current_plot_type|name, type='response')
pairs(est)

#### now i'll do a by-hand hurdle model on my own using a truncated negative binomial

### zeros and ones 
countdat$presence<-ifelse(countdat$value==0, 0, 1)
zerofit<-glmmTMB(presence~name*current_plot_type+(1 | block), family=binomial, data=countdat, REML=FALSE)
emmip(zerofit,~current_plot_type|name, type='response',CI=T)
est<-emmeans(zerofit, ~current_plot_type|name, type='response')
pairs(est)

### abundance with a truncated negbinom 
countdat$posicounts<-as.numeric(ifelse(countdat$value==0, "NA", countdat$value))
countfit<-glmmTMB(posicounts~name*current_plot_type+(1 | block), family=truncated_nbinom2(), data=countdat, REML=FALSE)
emmip(countfit,~current_plot_type|name, type='response',CI=T)
est<-emmeans(countfit, ~current_plot_type|name, type='response')
pairs(est)

############# treatment response: do analysis for total biomass by species #############

dat3<-(dat1[c(1:4,10,14,18,22)])
totwtdat<-as.data.frame(dat3 %>% pivot_longer(c(wt_max15_goro, wt_max15_tror, wt_max15_trcy)))
totwtdat$log_wt<-log(totwtdat$value) # log transform the weight data to get it normal looking

# model
totwtmod<-lmer(log_wt~name*current_plot_type+(1|block), data=totwtdat, REML=FALSE)

# test for fit, looks pretty good
sim<-simulateResiduals(totwtmod)
plot(sim)

# model summary
summary(totwtmod)
emmip(totwtmod,~current_plot_type|name, CI=T)

est<-emmeans(totwtmod,~current_plot_type|name, type='response')
pairs(est)

############# treatment response: do analysis for per capita biomass by species #############

dat4<-(dat1[c(1:4,9,13,17,22)])
pcwtdat<-as.data.frame(dat4 %>% pivot_longer(c(wt_percapita_goro, wt_percapita_tror, wt_percapita_trcy)))
pcwtdat$log_wt<-log(pcwtdat$value) # log transform the weight data to get it normal looking

# model
pcwtmod<-lmer(log_wt~name*current_plot_type+(1|block), data=pcwtdat, REML=FALSE)

# test for fit, looks pretty good
# sim<-simulateResiduals(pcwtmod)
# plot(sim)

# model summary
summary(pcwtmod)
emmip(pcwtmod,~current_plot_type|name,CI=T)

est<-emmeans(totwtmod,~current_plot_type|name, type='response')
pairs(est)

######### What about the log legacy - initial treatment for 2021 #######

############# log legacy response: do analysis for count by species ############# 
fit3_leg<-glmmTMB(value~name*initial+(1 | block), ziformula=~., family=nbinom2(), data=countdat, REML=FALSE)

## This stuff will just give you the end result counts with the zeros factored in...
summary(fit3_leg)
emmip(fit3_leg,~initial|name, type='response',CI=T)

est<-emmeans(fit3_leg,~initial|name, type='response')
pairs(est)

#### split up occurrence and abundance
### zeros and ones 
zerofit_leg<-glmmTMB(presence~name*initial+(1 | block), family=binomial, data=countdat, REML=FALSE)
emmip(zerofit_leg,~initial|name, type='response',CI=T)
est<-emmeans(zerofit_leg, ~initial|name, type='response')
pairs(est)

### abundance with a truncated negbinom 
countfit_leg<-glmmTMB(posicounts~name*initial+(1 | block), family=truncated_nbinom2(), data=countdat, REML=FALSE)
emmip(countfit_leg,~initial|name, type='response',CI=T)
est<-emmeans(countfit_leg, ~initial|name, type='response')
pairs(est)


############# log legacy response: do analysis for total weight by species ############# 

# weights
totwtmod_leg<-lmer(log_wt~name*initial+(1|block), data=totwtdat, REML=FALSE)

# model summary
summary(totwtmod_leg)
emmip(totwtmod_leg,~initial|name,CI=T)

est<-emmeans(totwtmod_leg,~initial|name, type='response')
pairs(est)

############# log legacy response: do analysis for per capita weight by species ############# 

# model
pcwtmod_leg<-lmer(log_wt~name*initial+(1|block), data=pcwtdat, REML=FALSE)

# test for fit, looks pretty good
# sim<-simulateResiduals(pcwtmod_leg)
# plot(sim)

# model summary
summary(pcwtmod_leg)
emmip(pcwtmod_leg,~initial|name,CI=T)

est<-emmeans(pcwtmod_leg,~initial|name, type='response')
pairs(est)


######### What about the physical barrier treatment for 2021 #######

############# physical barrier response: do analysis for count by species ############# 
fit3_phys<-glmmTMB(value~name*physical_barrier+(1 | block), ziformula=~., family=nbinom2(), data=countdat, REML=FALSE)

## This stuff will just give you the end result counts with the zeros factored in...
summary(fit3_phys)
emmip(fit3_phys,~physical_barrier|name, type='response',CI=T)

est<-emmeans(fit3_phys,~physical_barrier|name, type='response')
pairs(est)

#### split up occurrence and abundance
### zeros and ones 
zerofit_phys<-glmmTMB(presence~name*physical_barrier+(1 | block), family=binomial, data=countdat, REML=FALSE)
emmip(zerofit_phys,~physical_barrier|name, type='response',CI=T)
est<-emmeans(zerofit_phys, ~physical_barrier|name, type='response')
pairs(est)

### abundance with a truncated negbinom 
countfit_phys<-glmmTMB(posicounts~name*physical_barrier+(1 | block), family=truncated_nbinom2(), data=countdat, REML=FALSE)
emmip(countfit_phys,~physical_barrier|name, type='response',CI=T)
est<-emmeans(countfit_phys, ~physical_barrier|name, type='response')
pairs(est)

############# physical barrier response: do analysis for total weight by species ############# 

# weights
totwtmod_phys<-lmer(log_wt~name*physical_barrier+(1|block), data=totwtdat, REML=FALSE)

# model summary
summary(totwtmod_phys)
emmip(totwtmod_phys,~physical_barrier|name,CI=T)

est<-emmeans(totwtmod_phys,~physical_barrier|name, type='response')
pairs(est)

############# physical barrier response: do analysis for per capita weight by species ############# 

# model
pcwtmod_phys<-lmer(log_wt~name*physical_barrier+(1|block), data=pcwtdat, REML=FALSE)

# test for fit, looks pretty good
# sim<-simulateResiduals(pcwtmod_leg)
# plot(sim)

# model summary
summary(pcwtmod_phys)
emmip(pcwtmod_phys,~physical_barrier|name,CI=T)

est<-emmeans(pcwtmod_phys,~physical_barrier|name, type='response')
pairs(est)

#### What about both physical and legacy? #### 
### zeros and ones 
zerofit_intxn<-glmmTMB(presence~name*physical_barrier*initial+(1 | block), family=binomial, data=countdat, REML=FALSE)

sim<-simulateResiduals(zerofit_intxn)
plot(sim)
testDispersion(sim)
testZeroInflation(sim)

emmip(zerofit_intxn,initial~physical_barrier|name, type='response',CI=T)
est<-emmeans(zerofit_intxn, ~physical_barrier|name|initial, type='response')
pairs(est)

### abundance with a truncated negbinom 
countfit_intxn<-glmmTMB(posicounts~name*physical_barrier*initial+(1 | block), family=truncated_nbinom2(), data=countdat, REML=FALSE)

sim<-simulateResiduals(countfit_intxn)
plot(sim)
testDispersion(sim)
testZeroInflation(sim)

emmip(countfit_intxn,~physical_barrier|name|initial, type='response',CI=T)
est<-emmeans(countfit_intxn, ~physical_barrier|name|initial, type='response')
pairs(est)

### additive examp,es
zerofit_add<-glmmTMB(presence~name*physical_barrier+name*initial+(1 | block), family=binomial, data=countdat, REML=FALSE)

countfit_add<-glmmTMB(posicounts~name*physical_barrier+name*initial+(1 | block), family=truncated_nbinom2(), data=countdat, REML=FALSE)

### per capita biomass 

# model - intxn
pcwtmod_intxn<-lmer(log_wt~name*physical_barrier*initial+(1|block), data=pcwtdat, REML=FALSE)

# model - no intxn
pcwtmod_add<-lmer(log_wt~name*physical_barrier+name*initial+(1|block), data=pcwtdat, REML=FALSE)


#### 2021 model comparison for counts #### 
## zeros 
zero_candmods<-list("Plot type"=zerofit, 
                    "Physical barrier"=zerofit_phys,
                    "Nutrient island"=zerofit_leg,
                    "Physical Barrier + Nutrient Island"=zerofit_add,
                    "Physical Barrier x Nutrient Island"=zerofit_intxn)
aictab(zero_candmods)
pairs(emmeans(zerofit_add, ~physical_barrier|name|initial))


# counts
count_candmods<-list("Plot type"=countfit, 
                     "Physical barrier"=countfit_phys,
                     "Nutrient island"=countfit_leg,
                     "Physical Barrier + Nutrient Island"=countfit_add,
                     "Physical Barrier x Nutrient Island"=countfit_intxn)
aictab(count_candmods)
pairs(emmeans(countfit_leg, ~initial|name))


#### 2021 model comparison for per capita biomass #### 
pcwt_candmods<-list("Plot type"=pcwtmod, 
                    "Physical barrier"=pcwtmod_phys,
                    "Nutrient island"=pcwtmod_leg,
                    "Physical Barrier + Nutrient Island" = pcwtmod_add,
                    "Physical Barrier x Nutrient Island"=pcwtmod_intxn)
aictab(pcwt_candmods)
pairs(emmeans(pcwtmod_phys, ~physical_barrier|name))



### figures ###
# best fit models - zeros 
# colors
mimiscols<-c("#D66972","#108780")
# zerofit_add<-glmmTMB(presence~name*physical_barrier+name*initial+(1 | block), family=binomial, data=countdat, REML=FALSE)

zfit_est<-as.data.frame(emmeans(zerofit_add,~initial|physical_barrier|name, type='response'))
zfit_est$name<-c(rep("ngoro", 4), rep("ntrcy",4), rep("ntror",4))

pl1<-ggplot(zfit_est,aes(physical_barrier, prob,group=initial),)+
  scale_color_manual(values=mimiscols)+
  geom_point(aes(col=initial), size=2, position=position_dodge(width=0.5))+
  geom_linerange(aes(ymin=lower.CL, ymax=upper.CL, col=initial), position=position_dodge(width=0.5))+
  geom_line(aes(col=initial), position=position_dodge(width=0.5))+
  facet_wrap(vars(name))+
  theme_bw()+
  theme(strip.text.x = element_text(size=0),
        strip.background = element_blank(),
        axis.text=element_text(size=15),
        axis.title=element_text(size=20),
        legend.text=element_text(size=15),
        legend.title=element_text(size=15),
        legend.position="top")+
  xlab("Physical Barrier")+
  ylab("Occurrence")+
  geom_jitter(data=countdat,
             aes(x=physical_barrier, y=presence, color=initial), 
             height=0.1,
             alpha=0.5)+
  labs(color = "Initial Plot Type")
  
# best fit model - counts 
# countfit_leg<-glmmTMB(posicounts~name*initial+(1 | block), family=truncated_nbinom2(), data=countdat, REML=FALSE)

cfit_est<-as.data.frame(emmeans(countfit_leg,~initial|name, type='response'))
cfit_est$name<-c(rep("ngoro", 2), rep("ntrcy", 2), rep("ntror",2))

pl2<-ggplot(cfit_est,aes(initial, response, group=1))+
  geom_jitter(data=countdat,
              aes(x=initial, y=posicounts),
              width=0.1,
              alpha=0.4,
              color="gray")+
  geom_point(size=2, color="black")+
  geom_linerange(aes(ymin=lower.CL, ymax=upper.CL),color="black")+
  geom_line(color="black")+
  facet_wrap(vars(name))+
  theme_bw()+
  theme(strip.text.x = element_text(size=0),
        strip.background = element_blank(),
        axis.text=element_text(size=15),
        axis.title=element_text(size=20),
        legend.text=element_text(size=15),
        legend.title=element_text(size=15))+
  xlab("Initial Plot Type")+
  ylab("Abundance")+
  labs(color = "Initial Plot \n Type")

# best fit model - pcbiomass
# pcwtmod_phys<-lmer(log_wt~name*physical_barrier+(1|block), data=pcwtdat, REML=FALSE)

pcbfit_est<-as.data.frame(emmeans(pcwtmod_phys,~physical_barrier|name, type='response'))
#pcbfit_est$name<-c(rep("ngoro", 2), rep("ntrcy", 2), rep("ntror",2))


pl3<-ggplot(pcbfit_est,aes(physical_barrier, emmean, group=1))+
  geom_jitter(data=pcwtdat,
              aes(x=physical_barrier, y=log_wt),
              width=0.1,
              alpha=0.4,
              color="gray")+
  geom_point(size=2, color="black")+
  geom_linerange(aes(ymin=lower.CL, ymax=upper.CL),color="black")+
  facet_wrap(vars(name))+
  geom_line(color="black")+
  theme_bw()+
  theme(strip.text.x = element_text(size=0),
        strip.background = element_blank(),
        axis.text=element_text(size=15),
        axis.title=element_text(size=20),
        legend.text=element_text(size=15),
        legend.title=element_text(size=15))+
  xlab("Physical Barrier")+
  ylab("log(Biomass)")
pl3
### altogether now
ggarrange(pl1, pl2, pl3, ncol=1, common.legend = T)
pl1
pl2
pl3
