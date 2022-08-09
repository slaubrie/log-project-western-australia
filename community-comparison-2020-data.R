### code to make community figure 

#### packages #### 
require(vegan)
require(dplyr)
require(tidyr)
require(labdsv)
require(stringr)
require(ggplot2)
require(ggrepel)
require(lme4)
require(emmeans)


#### wrangling ####  
comm<-read.csv("species_composition_data_no_unk.csv", header=T)

# remove the locations surveyed in 2021 that were not surveyed in 2020 - aka, cm=0 and cm=21
comm<-comm[which(comm$cm_location!=21 & comm$cm_location!=0),]

# make a group name for each row
comm$grp<-apply(comm[c(1,3,5,6)], 1, paste, collapse=":")

# need to make each row a community using matrify
commsub<-comm[,c(15,10,13)] # group, species_code, and count of each species for each transect. transects are rows.
commtry<-matrify(commsub) # make it an expanded species matrix 
commtry$x[which(commtry$x=="x")]<-1 # "x" means that there were no individuals in the transect, but we are going to keep track of this as if it were a species
ncol(commtry) # how many species are we working with in our community matrix

# store grouping row names as a column, then remove rownames.
commtry$grps<-rownames(commtry)
rownames(commtry)<-NULL
names(commtry)

# split group info into columns for each variable
mat<-separate(commtry, 66, c("time","block","transect","init"), ":")
names(mat) #check

# add groupname using time, block, init columns
mat$grp<-apply(mat[c(66:69)], 1, paste, collapse=":")
names(mat) #check

# another df where the grouping variables are time, block, transect and initial state
# each row is a transect in a certain year.
df<-mat[,c(1:65,70)]
df2 = df %>% mutate(across(.cols=1:65,.fns=as.numeric)) # make everything numeric
rownames(df2)<-NULL # remove rownames

## new with group vars
nublock<-separate(df2, 66, c("time","block","transect", "init"), ":")

# want to sum across transects in same block X init X time treatment
nublock$sumgrp<-apply(mat[c(66,67,69)], 1, paste, collapse=":")
head(nublock)

# sum observations across initial X time X  block (group variable)
# this gives number of plants in each transect TYPE for each year in each block. should be 2 types X 2 years X 7 blocks rows 
blocksum<-rowsum(nublock[,c(1:65)], group=nublock$sumgrp)
blocksum$grps<-rownames(blocksum)
rownames(blocksum)<-NULL # remove rownames
nrow(blocksum) # it is 28 rows as expected

##  expand again
blocksum<-separate(blocksum, 66, c("time","block", "init"), ":")

# count info - just look at assemblies at initial timepoint.
# at the moment this includes where there were no plants ("x" column in matrix)
assemblies_t0<-blocksum[which(blocksum$time=="t0"),c(1:65)]

# group - these are the treatment variables that need to be separately fed into the MDS analaysis from the community analysis.
group_init<-blocksum$init[which(blocksum$time=='t0')]
group_block<-blocksum$block[which(blocksum$time=='t0')]

# MDS 
ass.rel.t0<-decostand(assemblies_t0, method='hel') #standardize assemblies 
ass.rel.t0_NMS<-metaMDS(ass.rel.t0, distance='bray', k=5, zerodist='ignore') # run MDS 
stressplot(ass.rel.t0_NMS) # check fit

# scores
mds_scores<-as.data.frame(scores(ass.rel.t0_NMS)$sites) # extract scores
mds_scores$site<-rownames(scores(ass.rel.t0_NMS)$sites) # extract names 
mds_scores$treatment<-group_init # grouping factor 1 
mds_scores$block<-group_block # grouping factor 2 

# explaining factors
init<-as.factor(group_init) # grouping factor 1- convert to factor
block<-as.factor(group_block) # grouping factor 2- convert to factor

### extracting species scores and plotting 
# species scores
species.scores<-as.data.frame(scores(ass.rel.t0_NMS,"species")) ## some species don't have scores
species.scores$species<-rownames(species.scores) 

### NMDS 1 and 2 
log<-mds_scores[mds_scores$treatment == "log", ][chull(mds_scores[mds_scores$treatment == 
                                                                    "log", c("NMDS1", "NMDS2")]), ]

open<-mds_scores[mds_scores$treatment == "open", ][chull(mds_scores[mds_scores$treatment == 
                                                                      "open", c("NMDS1", "NMDS2")]), ]

hulldat<-rbind(log,open)

ggplot()+
  theme_bw()+
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank(), 
        axis.text = element_text(size = 15),
        axis.title=element_text(size=20),
        legend.title=element_text(size=20), 
        legend.text=element_text(size=15))+
  geom_text_repel(data=species.scores, aes(NMDS1, NMDS2, label=species), alpha=0.9, size=5, col='darkgray', 
                  max.overlaps = getOption("ggrepel.max.overlaps", default = 30))+
  geom_polygon(data=hulldat, aes(NMDS1, NMDS2, fill=treatment, group=treatment), alpha=0.3)+scale_fill_manual(values=c("#D66972","#108780"), name="Initial \n Plot Type")+
  geom_point(data=mds_scores, aes(NMDS1, NMDS2, shape=block, col=treatment), size=6)+ scale_shape_manual(values = c(14,15,16,17,11,18,8), name='Block')+
  scale_colour_manual(values=c("#D66972","#108780"), name="Initial \n Plot Type")
  


# ### NMDS 3 and 4 
# log2<-mds_scores[mds_scores$treatment == "log", ][chull(mds_scores[mds_scores$treatment == 
#                                                                     "log", c("NMDS3", "NMDS4")]), ]
# 
# open2<-mds_scores[mds_scores$treatment == "open", ][chull(mds_scores[mds_scores$treatment == 
#                                                                       "open", c("NMDS3", "NMDS4")]), ]
# 
# hulldat2<-rbind(log2,open2)
# 
# 
# ggplot()+
#   theme_bw()+
#   theme(panel.background = element_blank(),
#         panel.grid.major = element_blank(),  #remove major-grid labels
#         panel.grid.minor = element_blank(),  #remove minor-grid labels
#         plot.background = element_blank(), 
#         axis.text = element_text(size = 12),
#         axis.title=element_text(size=15))+
#   geom_text_repel(data=species.scores, aes(NMDS3, NMDS4, label=species), alpha=0.9, size=3, col='darkgray')+
#   geom_polygon(data=hulldat2, aes(NMDS3, NMDS4, fill=treatment, group=treatment), alpha=0.3)+scale_fill_manual(values=c("#63A088","#56638A"), name="Treatment")+
#   geom_point(data=mds_scores, aes(NMDS3, NMDS4, shape=block, col=treatment), size=4)+ scale_shape_manual(values = c(14,15,16,17,11,18,8), name='Block')+
#   scale_colour_manual(values=c("#63A088","#56638A"), name="Treatment")



### varpart analysis
# can model using varpart to look at contributions of initial treatment and block
# block explain variation in community
trt_tot<-rda(ass.rel.t0~init+block) # run model using standardized data 
anova.cca(trt_tot) ## test for model significance
var.mod<-varpart(ass.rel.t0, init,block) # run model on standardized data
showvarparts(2, bg = c("#206713","#013F51"))

plot(var.mod, bg=c("#206713","#013F51"))
mtext("X1=Initial Plot Type; X2=Block", side=3)
## can test for significance of contribution of the fraction of initial treatment
# do this with partial redundancy analysis
trt_Frac<-rda(ass.rel.t0, init, block) # partial rda model
anova.cca(trt_Frac) ## this tells us if first condition, init, significantly contributes to overall variance explanation. 


#### Abundance analysis



#### wrangling ####  
comm<-read.csv("species_composition_data.csv", header=T)

# remove the locations surveyed in 2021 that were not surveyed in 2020 - aka, cm=0 and cm=21
comm<-comm[which(comm$cm_location!=21 & comm$cm_location!=0),]

# make a group name for each row
comm$grp<-apply(comm[c(1,3,5,6)], 1, paste, collapse=":")

# need to make each row a community using matrify
commsub<-comm[,c(15,10,13)] # group, species_code, and count of each species for each transect. transects are rows.
commtry<-matrify(commsub) # make it an expanded species matrix 
commtry$x[which(commtry$x=="x")]<-1 # "x" means that there were no individuals in the transect, but we are going to keep track of this as if it were a species
ncol(commtry) # how many species are we working with in our community matrix

# store grouping row names as a column, then remove rownames.
commtry$grps<-rownames(commtry)
rownames(commtry)<-NULL
names(commtry)

# split group info into columns for each variable
mat<-separate(commtry, 81, c("time","block","transect","init"), ":")
names(mat) #check

# add groupname using time, block, init columns
mat$grp<-apply(mat[c(81:84)], 1, paste, collapse=":")
names(mat) #check

# another df where the grouping variables are time, block, transect and initial state
# each row is a transect in a certain year.
df<-mat[,c(1:80,85)]
df2 = df %>% mutate(across(.cols=1:80,.fns=as.numeric)) # make everything numeric
rownames(df2)<-NULL # remove rownames

# sum observations across initial X transect X time X  block (group variable)
# this gives number of plants in each row observation
blocksum<-rowsum(df2[,c(1:80)], group=df2$grp)
blocksum$grps<-rownames(blocksum)
rownames(blocksum)<-NULL # remove rownames

## add in group vars 
nublock<-separate(blocksum, 81, c("time","block","transect","init"), ":")
nublock$total<-rowSums(nublock[,c(1:79)])
nublock$presence<-ifelse(nublock$total > 0,  1, 0)

# only before treatments installed
dat<-nublock[which(nublock$time=="t0"),]


# look at plant abundance in log vs open 

# look at range of data - what family should i use? 
range(dat$total)

## using poisson. including random effects throws is.Singular error: see ?is.Singular for details
abun.mod<-glm(total~init, data=dat, family='poisson') 
summary(abun.mod)

# the difference between open and log environments is marginal, where log has marginally significantly more plants. it's not a lot though. The model-estimated difference is like half a plant. 
emmeans(abun.mod, ~init, type='response')

# here's a quick and dirty plot of model estimated means and CIs for abundance in each kind of plot.
emmip(abun.mod, ~init, type='response', CI=T)+theme_bw()+labs(x="Initial condition", y="Number of plants")

# can also look at presence absence but it's really only a few zeros - no significant difference here.
abun.mod2<-glm(presence~init, data=dat, family='binomial')
summary(abun.mod2)
