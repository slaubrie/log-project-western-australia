### doing vital rate analysis

# packages 
require(lme4)
require(emmeans)
require(pscl)
require(glmmTMB)
require(tidyr)
require(DHARMa)

############## 2021 data #############

# read csv
dat<-read.csv('nplants_data_2021_git.csv', header=T)
dat1<-dat[which(dat$seeding_trt==1),]
dat1$physical_barrier<-as.factor(dat1$physical_barrier)
dat1$block<-as.factor(dat1$block)

############# treatment response: do analysis for count by species ############# 

dat2<-(dat1[c(1:4,7,11,15,22)])
countdat<-as.data.frame(dat2 %>% pivot_longer(c(ntror, ngoro, ntrcy)))
countmod<-glmmTMB(value~name*current_plot_type+(1|block), family="poisson", data=countdat)

# test for fit and zero inflation
# sim<-simulateResiduals(countmod)
# testZeroInflation(sim)
# plot(sim)

## going to do a hurdle model, which assumes a zero is only generated in one way 
# https://jsdajournal.springeropen.com/articles/10.1186/s40488-021-00121-4
# #https://stats.stackexchange.com/questions/81457/what-is-the-difference-between-zero-inflated-and-hurdle-models

# hurdle model 
# using examples as presented here: https://www.biorxiv.org/content/biorxiv/suppl/2017/05/01/132753.DC1/132753-2.pdf 
fit3<-glmmTMB(value~name*current_plot_type+(1 | block), ziformula=~., family=truncated_nbinom2(), data=countdat)
# sim3<-simulateResiduals(fit3)
# plot(sim3)
# testDispersion(sim3)
# testZeroInflation(sim3)

## This stuff will just give you the end result counts with the zeros factored in...
summary(fit3)
emmip(fit3,name~current_plot_type, type='response',CI=T)

est<-emmeans(fit3,~name|current_plot_type, type='response')
pairs(est)


#### now i'll do a by-hand hurdle model on my own using a truncated negative binomial

### zeros and ones 
countdat$presence<-ifelse(countdat$value==0, 0, 1)
zerofit<-glmmTMB(presence~name*current_plot_type+(1 | block), family=binomial, data=countdat)
emmip(zerofit,~name|current_plot_type, type='response',CI=T)
emmip(zerofit,~current_plot_type|name, type='response',CI=T)
est<-emmeans(zerofit, ~current_plot_type|name, type='response')
pairs(est)

### abundance with a truncated negbinom 
countdat$posicounts<-as.numeric(ifelse(countdat$value==0, "NA", countdat$value))
countfit<-glmmTMB(posicounts~name*current_plot_type+(1 | block), family=truncated_nbinom2(), data=countdat)
emmip(countfit,~current_plot_type|name, type='response',CI=T)
est<-emmeans(countfit, ~current_plot_type|name, type='response')
pairs(est)

############# treatment response: do analysis for total biomass by species #############

dat3<-(dat1[c(1:4,10,14,18,22)])
totwtdat<-as.data.frame(dat3 %>% pivot_longer(c(wt_max15_goro, wt_max15_tror, wt_max15_trcy)))
totwtdat$log_wt<-log(totwtdat$value) # log transform the weight data to get it normal looking

# model
totwtmod<-lmer(log_wt~name*current_plot_type+(1|block), data=totwtdat)

# test for fit, looks pretty good
sim<-simulateResiduals(totwtmod)
plot(sim)

# model summary
summary(totwtmod)
emmip(totwtmod,~current_plot_type|name, CI=T)

est<-emmeans(totwtmod,~name|current_plot_type, type='response')
pairs(est)

############# treatment response: do analysis for per capita biomass by species #############

dat4<-(dat1[c(1:4,9,13,17,22)])
pcwtdat<-as.data.frame(dat4 %>% pivot_longer(c(wt_percapita_goro, wt_percapita_tror, wt_percapita_trcy)))
pcwtdat$log_wt<-log(pcwtdat$value) # log transform the weight data to get it normal looking

# model
pcwtmod<-lmer(log_wt~name*current_plot_type+(1|block), data=pcwtdat)

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
fit3_leg<-glmmTMB(value~name*initial+(1 | block), ziformula=~., family=truncated_nbinom2(), data=countdat)

## This stuff will just give you the end result counts with the zeros factored in...
summary(fit3_leg)
emmip(fit3_leg,name~initial, type='response',CI=T)

est<-emmeans(fit3_leg,~initial|name, type='response')

#### split up occurrence and abundance
### zeros and ones 
zerofit_leg<-glmmTMB(presence~name*initial+(1 | block), family=binomial, data=countdat)
emmip(zerofit_leg,~name|initial, type='response',CI=T)
emmip(zerofit_leg,~initial|name, type='response',CI=T)
est<-emmeans(zerofit_leg, ~initial|name, type='response')
pairs(est)

### abundance with a truncated negbinom 
countfit_leg<-glmmTMB(posicounts~name*initial+(1 | block), family=truncated_nbinom2(), data=countdat)
emmip(countfit_leg,~initial|name, type='response',CI=T)
est<-emmeans(countfit_leg, ~initial|name, type='response')
pairs(est)

############# log legacy response: do analysis for total weight by species ############# 

# weights
totwtmod_leg<-lmer(log_wt~name*initial+(1|block), data=totwtdat)

# model summary
summary(totwtmod_leg)
emmip(totwtmod_leg,~initial|name,CI=T)

est<-emmeans(totwtmod_leg,~initial|name, type='response')
pairs(est)

############# log legacy response: do analysis for per capita weight by species ############# 

# model
pcwtmod_leg<-lmer(log_wt~name*initial+(1|block), data=pcwtdat)

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
fit3_phys<-glmmTMB(value~name*physical_barrier+(1 | block), ziformula=~., family=truncated_nbinom2(), data=countdat)

## This stuff will just give you the end result counts with the zeros factored in...
summary(fit3_phys)
emmip(fit3_phys,~physical_barrier|name, type='response',CI=T)

est<-emmeans(fit3_phys,~physical_barrier|name, type='response')
pairs(est)

#### split up occurrence and abundance
### zeros and ones 
zerofit_phys<-glmmTMB(presence~name*physical_barrier+(1 | block), family=binomial, data=countdat)
emmip(zerofit_phys,~physical_barrier|name, type='response',CI=T)
est<-emmeans(zerofit_phys, ~physical_barrier|name, type='response')
pairs(est)

### abundance with a truncated negbinom 
countfit_phys<-glmmTMB(posicounts~name*physical_barrier+(1 | block), family=truncated_nbinom2(), data=countdat)
emmip(countfit_phys,~physical_barrier|name, type='response',CI=T)
est<-emmeans(countfit_phys, ~physical_barrier|name, type='response')
pairs(est)

############# physical barrier response: do analysis for total weight by species ############# 

# weights
totwtmod_phys<-lmer(log_wt~name*physical_barrier+(1|block), data=totwtdat)

# model summary
summary(totwtmod_phys)
emmip(totwtmod_phys,~physical_barrier|name,CI=T)

est<-emmeans(totwtmod_phys,~physical_barrier|name, type='response')
pairs(est)

############# physical barrier response: do analysis for per capita weight by species ############# 

# model
pcwtmod_phys<-lmer(log_wt~name*physical_barrier+(1|block), data=pcwtdat)

# test for fit, looks pretty good
# sim<-simulateResiduals(pcwtmod_leg)
# plot(sim)

# model summary
summary(pcwtmod_phys)
emmip(pcwtmod_phys,~physical_barrier|name,CI=T)

est<-emmeans(pcwtmod_phys,~physical_barrier|name, type='response')
pairs(est)

############## 2022 data ##############

# read csv
dat<-read.csv('nplants_data_2022.csv', header=T)
dat1<-dat[which(dat$seeding_trt==1),]
dat1$physical_barrier<-as.factor(dat1$physical_barrier)
dat1$block<-as.factor(dat1$block)
dat2<-dat1[,c(1,3:4,7:9,13)]
dat22<-as.data.frame(dat2 %>% pivot_longer(c(ntrcy_germ, ngoro_germ, ntror_germ)))
dat22$attempts<-rep(15, nrow(dat22))
dat22$fails<-dat22$attempts-dat22$value
dat22$presence<-ifelse(countdat$value==0, 0, 1)
dat22$posicounts<-as.numeric(ifelse(dat22$value==0, "NA", dat22$value))

############# treatment response: do analysis for germination by species ############# 

fit<-glmmTMB(value~name*current_plot_type+(1 | block), ziformula=~., family=truncated_poisson(), data=dat22)
 sim<-simulateResiduals(fit)
 plot(sim)
 testDispersion(sim)
 testZeroInflation(sim)
 
 summary(fit)
 emmip(fit, ~current_plot_type|name, type='response', CI=T)
 est<-emmeans(fit,~current_plot_type|name, type='response')
 pairs(est)
 
######### What about the log legacy - initial treatment for 2022 #######

fit_leg<-glmmTMB(value~name*initial+(1 | block), ziformula=~., family=truncated_poisson(), data=dat22)
sim<-simulateResiduals(fit_leg)
plot(sim)
testDispersion(sim)
testZeroInflation(sim)

summary(fit_leg)
emmip(fit_leg, ~initial|name, type='response', CI=T)
est<-emmeans(fit_leg,~initial|name, type='response')
pairs(est)

# split up occurrence and abundance
### zeros and ones 
zerofit_leg<-glmmTMB(presence~name*initial+(1 | block), family=binomial, data=dat22)
emmip(zerofit_leg,~initial|name, type='response',CI=T)
est<-emmeans(zerofit_leg, ~initial|name, type='response')
pairs(est)

### abundance with a truncated negbinom 
countfit_leg<-glmmTMB(posicounts~name*initial+(1 | block), family=truncated_nbinom2(), data=dat22)
emmip(countfit_leg,~initial|name, type='response',CI=T)
est<-emmeans(countfit_leg, ~initial|name, type='response')
pairs(est)

######### What about the physical barrier - initial treatment for 2022 #######

fit_phys<-glmmTMB(value~name*physical_barrier+(1 | block), ziformula=~., family=truncated_poisson(), data=dat22)
sim<-simulateResiduals(fit_phys)
plot(sim)
testDispersion(sim)
testZeroInflation(sim)

summary(fit_phys)
emmip(fit_phys, ~physical_barrier|name, type='response', CI=T)
est<-emmeans(fit_phys,~physical_barrier|name, type='response')
pairs(est)

# split up occurrence and abundance
### zeros and ones 
zerofit_phys<-glmmTMB(presence~name*physical_barrier+(1 | block), family=binomial, data=dat22)
emmip(zerofit_phys,~physical_barrier|name, type='response',CI=T)
est<-emmeans(zerofit_phys, ~physical_barrier|name, type='response')
pairs(est)

### abundance with a truncated negbinom 
countfit_phys<-glmmTMB(posicounts~name*physical_barrier+(1 | block), family=truncated_nbinom2(), data=dat22)
emmip(countfit_phys,~physical_barrier|name, type='response',CI=T)
est<-emmeans(countfit_phys, ~physical_barrier|name, type='response')
pairs(est)