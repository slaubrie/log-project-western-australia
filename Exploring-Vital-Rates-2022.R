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

############## 2022 data ##############

# read csv
dat<-read.csv('nplants_data_2022.csv', header=T)
dat1<-dat[which(dat$seeding_trt==1),]
dat1$physical_barrier<-as.factor(dat1$physical_barrier)
dat1$block<-as.factor(dat1$block)
names(dat1)

# subset
dat2<-dat1[,c(1,2:4,10:12,19)] # these are block, transect, initial, current_plot_type, ngoro_plants, ntrcy_plants, ntror_plants, physical_barrier
head(dat2)

# pivot
dat22<-as.data.frame(dat2 %>% pivot_longer(c(ntrcy_plants, ngoro_plants, ntror_plants)))
range(dat22$value)
dat22$value>15 # one sample is larger than 15, it is a tror.

# max out at 15 
dat22$value<-as.numeric(ifelse(dat22$value>15, 15, dat22$value))

#################################### TREATMENT RESPONSE PRELIM ANALYSIS ####################################
############# treatment response: do analysis for count by species ############# 

# I am going to do the analysis as in 2021 now, with final counts and biomass. 
countmod<-glmmTMB(value~name*current_plot_type+(1|block), family="poisson", data=dat22)

# test for fit and zero inflation
sim<-simulateResiduals(countmod)
testZeroInflation(sim) # zero-inflated so need something else! 
plot(sim) 

## going to do a hurdle model, which assumes a zero is only generated in one way 
# https://jsdajournal.springeropen.com/articles/10.1186/s40488-021-00121-4
# #https://stats.stackexchange.com/questions/81457/what-is-the-difference-between-zero-inflated-and-hurdle-models

# hurdle model 
# using examples as presented here: https://www.biorxiv.org/content/biorxiv/suppl/2017/05/01/132753.DC1/132753-2.pdf 
fit3<-glmmTMB(value~name+current_plot_type+(1|block), ziformula=~., family=nbinom2(), data=dat22) # model has non-positive definite hessian matrix when the model expression contains an interaction, so made it additive.
 sim3<-simulateResiduals(fit3)
 plot(sim3)
 testDispersion(sim3) # looks ok ! 
 testZeroInflation(sim3) # looks ok ! 

summary(fit3)
emmip(fit3,~current_plot_type|name, type='response',CI=T)

# visualize
est<-emmeans(fit3,~current_plot_type|name, type='response')
pairs(est)

#### now i'll do a by-hand hurdle model on my own using a truncated negative binomial

### zeros and ones 
dat22$presence<-ifelse(dat22$value==0, 0, 1)
zerofit<-glmmTMB(presence~name*current_plot_type+(1 | block), family=binomial, data=dat22, REML=FALSE)
emmip(zerofit,~current_plot_type|name, type='response',CI=T)
est<-emmeans(zerofit, ~current_plot_type|name, type='response')
pairs(est)

### abundance with a truncated negbinom 
dat22$posicounts<-as.numeric(ifelse(dat22$value==0, "NA", dat22$value))
countfit<-glmmTMB(posicounts~name*current_plot_type+(1 | block), family=truncated_nbinom2(), data=dat22, REML=FALSE)
emmip(countfit,~current_plot_type|name, type='response',CI=T)
est<-emmeans(countfit, ~current_plot_type|name, type='response')
pairs(est)

############# treatment response: do analysis for total biomass by species #############

dat3<-(dat1[c(1:4,19,21,23,25)])
totwtdat<-as.data.frame(dat3 %>% pivot_longer(c(wt_max15_goro, wt_max15_tror, wt_max15_trcy)))
totwtdat$log_wt<-log1p(totwtdat$value) # log transformation does not help the insane amount of heteroscedasticity here, neither does a square root transformation.

# model
totwtmod<-lm(log_wt~name*current_plot_type, data=totwtdat) # fit is singular when including a random effect for block so not doing that

# test for fit, looks pretty good
sim<-simulateResiduals(totwtmod) # insanely heteroscedastic, doesn't get better with sqrt transform
plot(sim)

# model summary - don't trust this the residual dispersion is fucked, commenting out.
# summary(totwtmod)
# emmip(totwtmod,~current_plot_type|name, CI=T)
# 
# est<-emmeans(totwtmod,~current_plot_type|name, type='response')
# pairs(est)

############# treatment response: do analysis for per capita biomass by species #############

# per capita might help with heteroscedasticity  
dat4<-(dat1[c(1:4,19,20,22,24)])
pcwtdat<-as.data.frame(dat4 %>% pivot_longer(c(wt_percapita_goro, wt_percapita_tror, wt_percapita_trcy)))
pcwtdat$log_wt<-log(pcwtdat$value) # log transform the weight data to get it normal looking - it is heteroscedastic otherwise

# model
pcwtmod<-lm(log_wt~name*current_plot_type, data=pcwtdat) # fit is singular with random effect.

# test for fit, looks pretty good
sim<-simulateResiduals(pcwtmod)
plot(sim) ## much better! 

# model summary
summary(pcwtmod)
emmip(pcwtmod,~current_plot_type|name,CI=T)

est<-emmeans(totwtmod,~current_plot_type|name, type='response')
pairs(est)

#################################### PHYSICAL BARRIER ANALYSIS ####################################

######### What about the physical barrier - initial treatment for 2022 #######

fit_phys<-glmmTMB(value~name*physical_barrier+(1 | block), ziformula=~., family=poisson(), data=dat22, REML=F) #
sim<-simulateResiduals(fit_phys)
plot(sim)
testDispersion(sim) # not quite overdispersed; also if i fit nbinom model doesn't converge 
testZeroInflation(sim)

summary(fit_phys)
emmip(fit_phys, ~physical_barrier|name, type='response', CI=T)
est<-emmeans(fit_phys,~physical_barrier|name, type='response')
pairs(est)

# split up occurrence and abundance
### zeros and ones 
zerofit_phys<-glmmTMB(presence~name*physical_barrier+(1 | block), family=binomial, data=dat22, REML=F)

sim<-simulateResiduals(zerofit_phys)
plot(sim)
testDispersion(sim)
testZeroInflation(sim)

emmip(zerofit_phys,~physical_barrier|name, type='response',CI=T)
est<-emmeans(zerofit_phys, ~physical_barrier|name, type='response')
pairs(est)

### abundance with a truncated negbinom 
countfit_phys<-glmmTMB(posicounts~name*physical_barrier+(1 | block), family=truncated_nbinom2(), data=dat22, REML=F)

sim<-simulateResiduals(countfit_phys)
plot(sim)
testDispersion(sim)
testZeroInflation(sim)

emmip(countfit_phys,~physical_barrier|name, type='response',CI=T)
est<-emmeans(countfit_phys, ~physical_barrier|name, type='response')
pairs(est)

#################################### LEGACY ANALYSIS ####################################


# COME BACK TO THIS
#################################### PHYSICAL BARRIER X LEGACY ANALYSIS ####################################
######### What about both physical and legacy? ######### 
### zeros and ones 
zerofit_intxn<-glmmTMB(presence~name*physical_barrier*initial+(1 | block), family=binomial, data=dat22, REML=F)

sim<-simulateResiduals(zerofit_intxn)
plot(sim)
testDispersion(sim)
testZeroInflation(sim)

emmip(zerofit_intxn,initial~physical_barrier|name, type='response',CI=T)
est<-emmeans(zerofit_intxn, ~physical_barrier|name|initial, type='response')
pairs(est)

### abundance with a truncated negbinom 
countfit_intxn<-glmmTMB(posicounts~name*physical_barrier*initial+(1 | block), family=truncated_nbinom2(), data=dat22, REML=F)

sim<-simulateResiduals(countfit_intxn)
plot(sim)
testDispersion(sim)
testZeroInflation(sim)

emmip(countfit_intxn,~physical_barrier|name|initial, type='response',CI=T)
est<-emmeans(countfit_intxn, ~physical_barrier|name|initial, type='response')
pairs(est)


### adding 
zerofit_add<-glmmTMB(presence~name*physical_barrier+name*initial+(1 | block), family=binomial, data=dat22, REML=F)

countfit_add<-glmmTMB(posicounts~name*physical_barrier+name*initial+(1 | block), family=truncated_nbinom2(), data=dat22, REML=F)


##### 2022 candidate model comparison for counts ##### 

## zeros 
zero_candmods<-list("Plot type"=zerofit, 
                    "Physical barrier"=zerofit_phys,
                    "Nutrient island"=zerofit_leg,
                    "Physical Barrier + Nutrient Island"=zerofit_add,
                    "Physical Barrier x Nutrient Island"=zerofit_intxn)
aictab(zero_candmods)

# counts
count_candmods<-list("Plot type"=countfit, 
                     "Physical barrier"=countfit_phys,
                     "Nutrient island"=countfit_leg,
                     "Physical Barrier + Nutrient Island"=countfit_add,
                     "Physical Barrier x Nutrient Island"=countfit_intxn)
aictab(count_candmods)

### figures ###
# best fit models - zeros 
# colors
mimiscols<-c("#D66972","#108780")
# zerofit_add<-glmmTMB(presence~name*physical_barrier+name*initial+(1 | block), family=binomial, data=countdat, REML=FALSE)

zfit_est<-as.data.frame(emmeans(zerofit_add,~initial|physical_barrier|name, type='response'))
zfit_est$name<-c(rep("ngoro_germ", 4), rep("ntrcy_germ",4), rep("ntror_germ",4))

pl4<-ggplot(zfit_est,aes(physical_barrier, prob,group=initial),)+
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
  geom_jitter(data=dat22,
              aes(x=physical_barrier, y=presence, color=initial), 
              height=0.1,
              alpha=0.5)+
  labs(color = "Initial Plot Type")
pl4
pairs(emmeans(zerofit_add, ~physical_barrier|name|initial))
# best fit model - counts 
# countfit_leg<-glmmTMB(posicounts~name*initial+(1 | block), family=truncated_nbinom2(), data=countdat, REML=FALSE)

cfit_est<-as.data.frame(emmeans(countfit_leg,~initial|name, type='response'))
cfit_est$name<-c(rep("ngoro_germ", 2), rep("ntrcy_germ", 2), rep("ntror_germ",2))

pl5<-ggplot(cfit_est,aes(initial, response, group=1))+
  geom_jitter(data=dat22,
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
pairs(emmeans(countfit_leg, ~initial|name))

pl5
