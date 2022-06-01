################################
# Hierarchical Bayesian Model ##
# Response Variable: Pres/Abs  #
# Scale: local                 #
# Size Class: all              #
################################

################################
# Phylogenetic Distinctiveness #
# MDNN                         #  
# MDNC                         #
################################

require(brms)
require(labeling)
require(ape)
require(coda)
require(dplyr)
#require(phytools)
require(ggplot2)
require(broom)
#require(here)
require(psych)
require(lme4)
require(car)
require(tidybayes)

#########
# Data ##
#########

load("Cleaned_Data/Local_All_RelAbund_PresAbs_NoMETPOL_Final.RData")

local_all_paa$PlotID<-as.factor(local_all_paa$PlotID)

local_all_paa$Plot_Area2<-as.factor(local_all_paa$Plot_Area)

colnames(local_all_paa)[17]<-"phylo"

local_all_paa$phylo<-as.character(local_all_paa$phylo)

local_all_paa$species<-local_all_paa$phylo

######################
# add in phylogeny   #
######################

tree<-local_tree_paa

#scale co-variance matrix

nspp<-length(unique(tree$tip.label))
Vphy <- vcv(tree)
Vphy <- Vphy/(det(Vphy)^(1/nspp))

###############
# MDWNC ########
###############

mdnc<-filter(local_all_paa,PDmeasure=="MDWNC")
mdnc$PDmeasure<-droplevels(mdnc$PDmeasure)

#############
# quick QC ##
#############

length(unique(mdnc$PlotID)) # 310
length(unique(mdnc$phylo)) #  41
length(unique(tree$tip.label)) #  41

# transform + z-transform variables

mdnc$MAPs<-scale(mdnc$MAP,center=TRUE,scale=TRUE)
mdnc$MATs<-scale(mdnc$MAT,center=TRUE, scale=TRUE)
mdnc$Elev_ms<-scale(mdnc$Elev_m,center=TRUE, scale=TRUE)

mdnc$AridInd<-mdnc$AridInd*-1
mdnc$AridInds<-scale(mdnc$AridInd,center=TRUE,scale=TRUE)

mdnc$VPDs<-scale(mdnc$VPD,center=TRUE,scale=TRUE)
mdnc$HIIs<-scale(mdnc$HII,center=TRUE,scale=TRUE)
mdnc$HFPs<-scale(mdnc$HFP,center=TRUE,scale=TRUE)
mdnc$SubstrateAge_range_lg<-log(mdnc$SubstrateAge_range)
mdnc$Soilz<-scale(mdnc$SubstrateAge_range_lg,center=TRUE,scale=TRUE)
# mdnc$SeedWeights<-scale(log(mdnc$SeedWeight_g),center=TRUE,scale=TRUE)
# mdnc$WDs<-scale(mdnc$meanWD,center=TRUE,scale=TRUE)
# mdnc$maxHTs<-scale(mdnc$maxHT,center=TRUE,scale=TRUE)
mdnc$PD.stands<-scale(mdnc$PD.stand,center=TRUE, scale=TRUE)

mdnc$WD.diff.lrrs<-scale(mdnc$WD.diff.lrr,center=TRUE, scale=TRUE)
mdnc$WD.wt.diff.lrrs<-scale(mdnc$WD.wt.diff.lrr,center=TRUE, scale=TRUE)
mdnc$maxDBH.diff.lrrs<-scale(mdnc$maxDBH.diff.lrr,center=TRUE, scale=TRUE)
mdnc$maxDBH.wt.diff.lrrs<-scale(mdnc$maxDBH.wt.diff.lrr,center=TRUE, scale=TRUE)
mdnc$SM.diff.lrrs<-scale(mdnc$SM.diff.lrr,center=TRUE, scale=TRUE)
mdnc$SM.wt.diff.lrrs<-scale(mdnc$SM.wt.diff.lrr,center=TRUE, scale=TRUE)

##############################
# corelations of predictors  #
##############################

mdnc<-ungroup(mdnc)
mdnc_c<-select(mdnc, MATs, AridInds,HIIs, Soilz,
               PD.stands, WD.wt.diff.lrrs,maxDBH.wt.diff.lrrs,
               SM.wt.diff.lrrs)

p_out<-corr.test(mdnc_c,method="pearson") # all below 0.7

####################
# model structure: #
####################

prior2 <- get_prior(PresAbs ~  AridInds+MATs+HIIs+Soilz+ PD.stands+WD.wt.diff.lrrs+maxDBH.wt.diff.lrrs+SM.wt.diff.lrrs+
                      PD.stands:MATs+ PD.stands:AridInds+PD.stands:Soilz+PD.stands:HIIs+
                      WD.wt.diff.lrrs:MATs+WD.wt.diff.lrrs:AridInds+WD.wt.diff.lrrs:Soilz+WD.wt.diff.lrrs:HIIs+
                      maxDBH.wt.diff.lrrs:MATs+maxDBH.wt.diff.lrrs:AridInds+maxDBH.wt.diff.lrrs:Soilz+maxDBH.wt.diff.lrrs:HIIs+ 
                      SM.wt.diff.lrrs:MATs+SM.wt.diff.lrrs:AridInds+SM.wt.diff.lrrs:Soilz+SM.wt.diff.lrrs:HIIs+ 
                      (1|Plot_Area2/PlotID)+ (1|species)+ (1|gr(phylo, cov=Vphy)),
                    data = mdnc, 
                    data2=list(Vphy=Vphy),
                    family = bernoulli())

# fit model

mdnc_bern2 <- brm(PresAbs ~  AridInds+MATs+HIIs+Soilz+ PD.stands+WD.wt.diff.lrrs+maxDBH.wt.diff.lrrs+SM.wt.diff.lrrs+
                    PD.stands:MATs+ PD.stands:AridInds+PD.stands:Soilz+PD.stands:HIIs+
                    WD.wt.diff.lrrs:MATs+WD.wt.diff.lrrs:AridInds+WD.wt.diff.lrrs:Soilz+WD.wt.diff.lrrs:HIIs+
                    maxDBH.wt.diff.lrrs:MATs+maxDBH.wt.diff.lrrs:AridInds+maxDBH.wt.diff.lrrs:Soilz+maxDBH.wt.diff.lrrs:HIIs+ 
                    SM.wt.diff.lrrs:MATs+SM.wt.diff.lrrs:AridInds+SM.wt.diff.lrrs:Soilz+SM.wt.diff.lrrs:HIIs+ 
                    (1|Plot_Area2/PlotID)+ (1|species)+ (1|gr(phylo, cov=Vphy)),
                  data = mdnc, 
                  data2=list(Vphy=Vphy),
                  family = bernoulli(), 
                  prior = prior2,
                  sample_prior = TRUE, chains = 4, cores =4, control=list(adapt_delta=0.99), 
                  iter = 4500, warmup = 2500)

######

pp_check(mdnc_bern2,nsamples=100) #+scale_x_continuous(trans="log")

np <- nuts_params(mdnc_bern2)
str(np)
# extract the number of divergence transitions
sum(subset(np, Parameter == "divergent__")$Value) #

summary(rhat(mdnc_bern2)) 

summary(mdnc_bern2,waic=FALSE, loo=FALSE)

plot(mdnc_bern2)

plot(conditional_effects(mdnc_bern2), points = FALSE) 

save(mdnc_bern2,file="~/MacroEcology_group/members/Dylan/Hawaii_Invasion/Cleaned_Data/BRMS_Local_All_PresAbs_MDNC_NoMETPOL_output.RData")

# remove non-significant interactions

bern_coef2<-mdnc_bern2 %>%
  #recover_types(rootsNSC_bm_simple$data)%>%
  spread_draws(b_Intercept, b_AridInds, b_MATs, b_HIIs, b_Soilz, b_PD.stands,b_WD.wt.diff.lrrs,
               b_maxDBH.wt.diff.lrrs,b_SM.wt.diff.lrrs,
               `b_MATs:PD.stands`,`b_AridInds:PD.stands`,`b_Soilz:PD.stands`,             
               `b_HIIs:PD.stands`,`b_MATs:WD.wt.diff.lrrs`,`b_AridInds:WD.wt.diff.lrrs`,    
               `b_Soilz:WD.wt.diff.lrrs`,`b_HIIs:WD.wt.diff.lrrs`,`b_MATs:maxDBH.wt.diff.lrrs`,    
               `b_AridInds:maxDBH.wt.diff.lrrs`,`b_Soilz:maxDBH.wt.diff.lrrs`,`b_HIIs:maxDBH.wt.diff.lrrs`,    
               `b_MATs:SM.wt.diff.lrrs`,`b_AridInds:SM.wt.diff.lrrs`,`b_Soilz:SM.wt.diff.lrrs`,       
               `b_HIIs:SM.wt.diff.lrrs`) %>%
  reshape2::melt(.,id.vars=c(".chain",".iteration",".draw"),variable.name="terms",na.rm=T)

bern_coef2<-dplyr::filter(bern_coef2, terms!="b_Intercept")

bern_coef2_pp<-ggplot(bern_coef2, aes(y=terms, x=value))+
  geom_vline(xintercept=0,color="gray60",linetype="solid",size=0.5) +
  geom_halfeyeh(point_interval = median_hdci, .width=c(0.68,0.95), 
                color="#2171b5", fill="#bdd7e7",size=0.7,alpha=0.5) +
  labs(y="",x="Median standardized estimate")+
  
  theme_bw()+theme(axis.title.y=element_text(colour="black",face="bold",size=8),
                   axis.title.x=element_text(colour="black",face="bold",size=8),
                   axis.text.x=element_text(colour="black",face="bold",size=8),
                   axis.text.y=element_text(colour=c("black"),face="bold",size=8),
                   plot.margin=margin(t=0.5,r=0.1,b=0.1,l=0.1,unit="cm"),
                   legend.position = "none",
                   panel.grid.major = element_blank(), panel.grid.minor = element_blank())


#significant interactions
# Soilz:SM.wt.diff.lrrs 
# MATs:maxDBH.wt.diff.lrrs 
# HIIs:WD.wt.diff.lrrs 

# then removed MATs:maxDBH.wt.diff.lrrs

# reduced model
prior3 <- get_prior(PresAbs ~ AridInds+MATs+HIIs+Soilz+ PD.stands+WD.wt.diff.lrrs+maxDBH.wt.diff.lrrs+SM.wt.diff.lrrs+
                      Soilz:SM.wt.diff.lrrs+HIIs:WD.wt.diff.lrrs +  
                      (1|Plot_Area2/PlotID)+ (1|species)+ (1|gr(phylo, cov=Vphy)),
                    data = mdnc, 
                    data2=list(Vphy=Vphy),
                    family = bernoulli())

# fit model

mdnc_bern2_red <- brm(PresAbs ~ AridInds+MATs+HIIs+Soilz+ PD.stands+WD.wt.diff.lrrs+maxDBH.wt.diff.lrrs+SM.wt.diff.lrrs+
                        Soilz:SM.wt.diff.lrrs+HIIs:WD.wt.diff.lrrs +  
                        (1|Plot_Area2/PlotID)+ (1|species)+ (1|gr(phylo, cov=Vphy)),
                      data = mdnc, 
                      data2=list(Vphy=Vphy),
                      family = bernoulli(), 
                      prior = prior3,
                      sample_prior = TRUE, chains = 4, cores =4, control=list(adapt_delta=0.9999, max_treedepth=12), 
                      iter = 4500, warmup = 2500)

pp_check(mdnc_bern2_red,nsamples=100) #+scale_x_continuous(trans="log")

np <- nuts_params(mdnc_bern2_red)
str(np)
# extract the number of divergence transitions
sum(subset(np, Parameter == "divergent__")$Value) #

summary(rhat(mdnc_bern2_red)) 

summary(mdnc_bern2_red,waic=FALSE, loo=FALSE)

# check 
bern_coef22<-mdnc_bern2_red %>%
  #recover_types(rootsNSC_bm_simple$data)%>%
  spread_draws(b_Intercept, b_AridInds, b_MATs, b_HIIs, b_Soilz, b_PD.stands,b_WD.wt.diff.lrrs,
               b_maxDBH.wt.diff.lrrs,b_SM.wt.diff.lrrs,
               `b_HIIs:WD.wt.diff.lrrs`,    
               `b_Soilz:SM.wt.diff.lrrs`) %>%
  reshape2::melt(.,id.vars=c(".chain",".iteration",".draw"),variable.name="terms",na.rm=T)

bern_coef22<-dplyr::filter(bern_coef22, terms!="b_Intercept")

bern_coef22_pp<-ggplot(bern_coef22, aes(y=terms, x=value))+
  geom_vline(xintercept=0,color="gray60",linetype="solid",size=0.5) +
  geom_halfeyeh(point_interval = median_hdci, .width=c(0.68,0.95), 
                color="#2171b5", fill="#bdd7e7",size=0.7,alpha=0.5) +
  labs(y="",x="Median standardized estimate")+
  
  theme_bw()+theme(axis.title.y=element_text(colour="black",face="bold",size=8),
                   axis.title.x=element_text(colour="black",face="bold",size=8),
                   axis.text.x=element_text(colour="black",face="bold",size=8),
                   axis.text.y=element_text(colour=c("black"),face="bold",size=8),
                   plot.margin=margin(t=0.5,r=0.1,b=0.1,l=0.1,unit="cm"),
                   legend.position = "none",
                   panel.grid.major = element_blank(), panel.grid.minor = element_blank())


save(mdnc_bern2, mdnc_bern2_red, file="~/MacroEcology_group/members/Dylan/Hawaii_Invasion/Cleaned_Data/BRMS_Local_All_PresAbs_MDNC_NoMETPOL_output.RData")
