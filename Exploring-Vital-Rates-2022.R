

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
# I think this doens't work because it's mostly zeros.
# fit<-glmmTMB(value~name*current_plot_type+(1| block), ziformula=~., family=poisson(), data=dat22)
#  sim<-simulateResiduals(fit)
#  plot(sim)
#  testDispersion(sim)
#  testZeroInflation(sim)
#  
#  summary(fit)
#  emmip(fit, ~current_plot_type|name, type='response', CI=T)
#  est<-emmeans(fit,~current_plot_type|name, type='response')
#  pairs(est)

### zeros and ones 
zerofit<-glmmTMB(presence~name*current_plot_type+(1 | block), family=binomial, data=dat22, REML=F)
emmip(zerofit,~current_plot_type|name, type='response',CI=T)
est<-emmeans(zerofit, ~current_plot_type|name, type='response')
pairs(est)

### abundance with a truncated negbinom 
countfit<-glmmTMB(posicounts~name*current_plot_type+(1 | block), family=truncated_nbinom2(), data=dat22, REML=F)
emmip(countfit,~current_plot_type|name, type='response',CI=T)
est<-emmeans(countfit, ~current_plot_type|name, type='response')
pairs(est)

######### What about the log legacy - initial treatment for 2022 #######

fit_leg<-glmmTMB(value~name*initial+(1 | block), ziformula=~., family=poisson(), data=dat22, REML=F)
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
zerofit_leg<-glmmTMB(presence~name*initial+(1 | block), family=binomial, data=dat22, REML=F)
emmip(zerofit_leg,~initial|name, type='response',CI=T)
est<-emmeans(zerofit_leg, ~initial|name, type='response')
pairs(est)

### abundance with a truncated negbinom 
countfit_leg<-glmmTMB(posicounts~name*initial+(1 | block), family=truncated_nbinom2(), data=dat22, REML=F)
emmip(countfit_leg,~initial|name, type='response',CI=T)+
  theme_bw()
est<-emmeans(countfit_leg, ~initial|name, type='response')
pairs(est)

######### What about the physical barrier - initial treatment for 2022 #######

fit_phys<-glmmTMB(value~name*physical_barrier+(1 | block), ziformula=~., family=nbinom2(), data=dat22, REML=F)
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


#### What about both physical and legacy? #### 
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
