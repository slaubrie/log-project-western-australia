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
dat2<-dat1[,c(1,2:4,26:28,19)] # these are block, transect, initial, current_plot_type, ngoro_plants, ntrcy_plants, ntror_plants, physical_barrier
head(dat2)

## because of the way weeding worked, where we thinned only once and probably too early in the season, there were cases where plants came up, we thinned them, and then when we collected plants at the end of the growing season there were more that popped up. 
# I am choosing to use the total number of a species that popped up in the zone of planting. this required a somewhat complex excel formula that did not double-count the one individual when summing the total, and also accounts for the plants that came up but
# did not survive. This value is called 'tot' (e.g. ngoro_tot). 
# It's impossible to get the same dataset as the one from 2021 because we planted the plants and only came back when it was time to collect at the end of the season, thus capturing both survival through the season and germination. 
# to get the total in the "tot" columns, I wrote a formula in excel, written below
# total = if((t0+t2=0),0,if(t0=0,t2,if(t2=0,t0, if(t0+t2=1,1,(t0-1+t2)))))
# in this expression t0 is the germ value (e.g. ngoro_germ) and t2 is the plants value (e.g. ngoro_plants). 
# I can run the analysis with germ, plants, or tot. I'm choosing to do the analysis on tot, but we can revisit this choice if we feel there's a better way to comparably (or not) analyze the data between 2021 and 2023. 
 
# pivot
dat22<-as.data.frame(dat2 %>% pivot_longer(c(ntrcy_tot, ngoro_tot, ntror_tot)))
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
# no significant differences

### abundance with a truncated negbinom 
dat22$posicounts<-as.numeric(ifelse(dat22$value==0, "NA", dat22$value))
countfit<-glmmTMB(posicounts~name*current_plot_type+(1 | block), family=truncated_nbinom2(), data=dat22, REML=FALSE)
emmip(countfit,~current_plot_type|name, type='response',CI=T)
est<-emmeans(countfit, ~current_plot_type|name, type='response')
pairs(est)
# no significant differences in abundance

############# treatment response: do analysis for total biomass by species #############


dat3<-(dat1[c(1:4,19,21,23,25)])
totwtdat<-as.data.frame(dat3 %>% pivot_longer(c(wt_max15_goro, wt_max15_tror, wt_max15_trcy)))
totwtdat$log_wt<-log1p(totwtdat$value) # log transformation does not help the insane amount of heteroscedasticity here, neither does a square root transformation.

# model
totwtmod<-lm(log_wt~name*current_plot_type, data=totwtdat) # fit is singular when including a random effect for block so not doing that

# test for fit
sim<-simulateResiduals(totwtmod) # insanely heteroscedastic, doesn't get better with sqrt transform (all transformations make it worse)
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
est2<-emmeans(totwtmod,~name, type='response')
pairs(est)
# not any differences in per capita biomass among treatments within species, but goro is larger than the other two

#################################### PHYSICAL BARRIER ANALYSIS ####################################

######### What about the physical barrier - initial treatment for 2022 #######

# this is a model where do zero inflation and count together (occurrence and abundance)
fit_phys<-glmmTMB(value~name*physical_barrier+(1 | block), ziformula=~., family=nbinom2(), data=dat22, REML=F) #
sim<-simulateResiduals(fit_phys)
plot(sim)
testDispersion(sim) 
testZeroInflation(sim)

summary(fit_phys)
emmip(fit_phys, ~physical_barrier|name, type='response', CI=T)
est<-emmeans(fit_phys,~physical_barrier|name, type='response')
pairs(est)
# no significant differences in total plants (occurrence and abundance) 

# split up occurrence and abundance
### zeros and ones 
zerofit_phys<-glmmTMB(presence~name*physical_barrier+(1 | block), family=binomial, data=dat22, REML=F)

sim<-simulateResiduals(zerofit_phys)
plot(sim)
testDispersion(sim) # looks good
testZeroInflation(sim) # looks good

summary(zerofit_phys)
emmip(zerofit_phys,~physical_barrier|name, type='response',CI=T)
est<-emmeans(zerofit_phys, ~physical_barrier|name, type='response')
pairs(est)

# physical barrier occurrence results
# physical barrier level does not significantly explain variation in goro 
# trcy and tror probability of occurrence are higher in places where there is a physical barrier (p=0.02 and p=0.009 respectively)

### abundance with a truncated negbinom 
countfit_phys<-glmmTMB(posicounts~name*physical_barrier+(1 | block), family=truncated_nbinom2(), data=dat22, REML=F)

sim<-simulateResiduals(countfit_phys)
plot(sim)
testDispersion(sim)
testZeroInflation(sim)

summary(countfit_phys)
emmip(countfit_phys,~physical_barrier|name, type='response',CI=T)
est<-emmeans(countfit_phys, ~physical_barrier|name, type='response')
pairs(est)
# abundance does not differ between physical barrier treatments for any of the plant species


# model
# pcwtmod_phys<-lmer(log_wt~name*physical_barrier+(1|block), data=pcwtdat, REML=FALSE) # singular 
pcwtmod_phys<-lm(log_wt~name*physical_barrier, data=pcwtdat) 


# test for fit, looks pretty good
 sim<-simulateResiduals(pcwtmod_phys)
 plot(sim)

# model summary
summary(pcwtmod_phys)
emmip(pcwtmod_phys,~physical_barrier|name,CI=T)

est<-emmeans(pcwtmod_phys,~physical_barrier|name, type='response')
pairs(est)

#################################### LEGACY ANALYSIS ####################################

############# log legacy response: do analysis for count by species ############# 
# fit3_leg<-glmmTMB(value~name*initial+(1 | block), ziformula=~., family=nbinom2(), data=dat22, REML=FALSE) # does not converge, i think it's poisson dist.
fit3_leg<-glmmTMB(value~name*initial+(1 | block), ziformula=~., family=poisson(), data=dat22, REML=FALSE) 

## This stuff will just give you the end result counts with the zeros factored in...
summary(fit3_leg)
emmip(fit3_leg,~initial|name, type='response',CI=T)

est<-emmeans(fit3_leg,~initial|name, type='response')
pairs(est)

#### split up occurrence and abundance
### zeros and ones 
zerofit_leg<-glmmTMB(presence~name*initial+(1 | block), family=binomial, data=dat22, REML=FALSE)
emmip(zerofit_leg,~initial|name, type='response',CI=T)
est<-emmeans(zerofit_leg, ~initial|name, type='response')
pairs(est)

# probability of occurrence is lower for trcy in places where initial treatment is open (higher where there is a log legacy p=0.03)

### abundance with a truncated negbinom 
countfit_leg<-glmmTMB(posicounts~name*initial+(1 | block), family=truncated_nbinom2(), data=dat22, REML=FALSE)
emmip(countfit_leg,~initial|name, type='response',CI=T)
est<-emmeans(countfit_leg, ~initial|name, type='response')
pairs(est)


# model for per capita weight
# pcwtmod_leg<-lmer(log_wt~name*initial+(1|block), data=pcwtdat, REML=FALSE) #singular
pcwtmod_leg<-lm(log_wt~name*initial, data=pcwtdat)

# test for fit, looks pretty good
# sim<-simulateResiduals(pcwtmod_leg)
# plot(sim)

# model summary
summary(pcwtmod_leg)
emmip(pcwtmod_leg,~initial|name,CI=T)

est<-emmeans(pcwtmod_leg,~initial|name, type='response')
pairs(est)

# biomass per capita in tror is lower in open legacy environments (p=0.03); biomass per capita in trcy is higher in open legacy environments

#################################### PHYSICAL BARRIER X LEGACY ANALYSIS ####################################
######### What about both physical and legacy? ######### 

#### interaction model
### zeros and ones 
zerofit_intxn<-glmmTMB(presence~name*physical_barrier*initial+(1 | block), family=binomial, data=dat22, REML=F)

sim<-simulateResiduals(zerofit_intxn)
plot(sim)
testDispersion(sim)
testZeroInflation(sim)

emmip(zerofit_intxn,initial~physical_barrier|name, type='response',CI=T)
est<-emmeans(zerofit_intxn, ~physical_barrier|name|initial, type='response')
pairs(est)

# goro presence/absence not explained by physical barrier or initial treatment 
# when the initial treatment is "open", trcy has higher prbability of occurring when there is a physical barrier. the same goes for tror.
# if there is a legacy of a log, then there is no significant difference between physical barrier treatments for tror or for trcy. 

### abundance with a truncated negbinom 
countfit_intxn<-glmmTMB(posicounts~name*physical_barrier*initial+(1 | block), family=truncated_nbinom2(), data=dat22, REML=F)

sim<-simulateResiduals(countfit_intxn)
plot(sim)
testDispersion(sim)
testZeroInflation(sim)

emmip(countfit_intxn,~physical_barrier|name|initial, type='response',CI=T)
est<-emmeans(countfit_intxn, ~physical_barrier|name|initial, type='response')
pairs(est)

# abundance not explained by physical barrier or initial treatment for all three species

#### additive model 
zerofit_add<-glmmTMB(presence~name*physical_barrier+name*initial+(1 | block), family=binomial, data=dat22, REML=F)
summary(zerofit_add)

countfit_add<-glmmTMB(posicounts~name*physical_barrier+name*initial+(1 | block), family=truncated_nbinom2(), data=dat22, REML=F)
summary(countfit_add)

### per capita biomass 

# model - intxn
# pcwtmod_intxn<-lmer(log_wt~name*physical_barrier*initial+(1|block), data=pcwtdat, REML=FALSE) # singular
pcwtmod_intxn<-lm(log_wt~name*physical_barrier*initial, data=pcwtdat)

# model - no intxn
# pcwtmod_add<-lmer(log_wt~name*physical_barrier+name*initial+(1|block), data=pcwtdat, REML=FALSE) # singular 
 pcwtmod_add<-lm(log_wt~name*physical_barrier+name*initial, data=pcwtdat) 
summary(pcwtmod_add)

emmip(pcwtmod_add,~physical_barrier|name|initial, type='response',CI=T)
est<-emmeans(pcwtmod_add, ~initial|name, type='response')
pairs(est)

##### 2022 candidate model comparison ##### 

## zeros 
zero_candmods<-list("Plot type"=zerofit, 
                    "Physical barrier"=zerofit_phys,
                    "Nutrient island"=zerofit_leg,
                    "Physical Barrier + Nutrient Island"=zerofit_add,
                    "Physical Barrier x Nutrient Island"=zerofit_intxn)
aictab(zero_candmods)

# best fit model is a tie between physical barrier and physical barrier + nutrient island

# counts
count_candmods<-list("Plot type"=countfit, 
                     "Physical barrier"=countfit_phys,
                     "Nutrient island"=countfit_leg,
                     "Physical Barrier + Nutrient Island"=countfit_add,
                     "Physical Barrier x Nutrient Island"=countfit_intxn)
aictab(count_candmods)

pcwt_candmods<-list("Plot type"=pcwtmod, 
                    "Physical barrier"=pcwtmod_phys,
                    "Nutrient island"=pcwtmod_leg,
                    "Physical Barrier + Nutrient Island" = pcwtmod_add,
                    "Physical Barrier x Nutrient Island"=pcwtmod_intxn)
aictab(pcwt_candmods)


# tror is bigger when the legacy of the log and the log is still present

### figures ###
# best fit models - zeros 
# colors
mimiscols<-c("#D66972","#108780")
zerofit_add<-glmmTMB(presence~name*physical_barrier+name*initial+(1 | block), family=binomial, data=dat22, REML=FALSE)

zfit_est<-as.data.frame(emmeans(zerofit_add,~initial|physical_barrier|name, type='response'))

pl4<-ggplot(zfit_est,aes(physical_barrier,prob,group=initial),)+
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
summary(zerofit_add)
pairs(emmeans(zerofit_add, ~physical_barrier|name|initial))

# tror and trcy do better where there is a physical barrier as compared to when there is not (tror p=0.03, trcy p=0.007)

pairs(emmeans(zerofit_add, ~initial|name|physical_barrier))

# trcy does better with log initial as compared to open initial (p=0.03)

# best fit model - counts 
 countfit_leg<-glmmTMB(posicounts~name*initial+(1 | block), family=truncated_nbinom2(), data=dat22, REML=FALSE)

cfit_est<-as.data.frame(emmeans(countfit_leg,~initial|name, type='response'))

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

# best fit model - pcbiomass

pcbfit_est<-as.data.frame(emmeans(pcwtmod_add,~physical_barrier|name|initial, type='response'))


pl6<-ggplot(pcbfit_est,aes(physical_barrier, emmean, group=initial))+
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
  ylab("log(Biomass)")+
  geom_jitter(data=pcwtdat,
              aes(x=physical_barrier, y=log_wt, color=initial), 
              height=0.1,
              alpha=0.5)+
  labs(color = "Initial Plot Type")
pl6
summary(pcwtmod_add)
pairs(emmeans(pcwtmod_add, ~initial|name))
pairs(emmeans(pcwtmod_add, ~physical_barrier|name))
