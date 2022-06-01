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
#require(geoR)
require(lme4)
require(car)
require(tidybayes)

#########
# Data ##
#########

load("~/MacroEcology_group/members/Dylan/Hawaii_Invasion/Cleaned_Data/Local_All_RelAbund_PresAbs_Final.RData")

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

length(unique(mdnc$PlotID)) # 460
length(unique(mdnc$phylo)) #  41
length(unique(tree$tip.label)) #  41

# transform + z-transform variables

mdnc$Island<-as.factor(mdnc$Island)

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

# model without plot area as random group
# get priors
# 

# model with plot area as random group

prior2 <- get_prior(PresAbs ~  AridInds+MATs+HIIs+Soilz+ PD.stands+WD.wt.diff.lrrs+maxDBH.wt.diff.lrrs+SM.wt.diff.lrrs+
                      PD.stands:MATs+ PD.stands:AridInds+PD.stands:Soilz+PD.stands:HIIs+
                      WD.wt.diff.lrrs:MATs+WD.wt.diff.lrrs:AridInds+WD.wt.diff.lrrs:Soilz+WD.wt.diff.lrrs:HIIs+
                      maxDBH.wt.diff.lrrs:MATs+maxDBH.wt.diff.lrrs:AridInds+maxDBH.wt.diff.lrrs:Soilz+maxDBH.wt.diff.lrrs:HIIs+ 
                      SM.wt.diff.lrrs:MATs+SM.wt.diff.lrrs:AridInds+SM.wt.diff.lrrs:Soilz+SM.wt.diff.lrrs:HIIs+ 
                      (1|Plot_Area2/PlotID)+ (1|Island)+(1|species)+ (1|gr(phylo, cov=Vphy)),
                    data = mdnc, 
                    data2=list(Vphy=Vphy),
                    family = bernoulli())

# fit model

mdnc_bern2 <- brm(PresAbs ~  AridInds+MATs+HIIs+Soilz+ PD.stands+WD.wt.diff.lrrs+maxDBH.wt.diff.lrrs+SM.wt.diff.lrrs+
                    PD.stands:MATs+ PD.stands:AridInds+PD.stands:Soilz+PD.stands:HIIs+
                    WD.wt.diff.lrrs:MATs+WD.wt.diff.lrrs:AridInds+WD.wt.diff.lrrs:Soilz+WD.wt.diff.lrrs:HIIs+
                    maxDBH.wt.diff.lrrs:MATs+maxDBH.wt.diff.lrrs:AridInds+maxDBH.wt.diff.lrrs:Soilz+maxDBH.wt.diff.lrrs:HIIs+ 
                    SM.wt.diff.lrrs:MATs+SM.wt.diff.lrrs:AridInds+SM.wt.diff.lrrs:Soilz+SM.wt.diff.lrrs:HIIs+ 
                    (1|Plot_Area2/PlotID)+ (1|Island)+(1|species)+ (1|gr(phylo, cov=Vphy)),
                    data = mdnc, 
                   data2=list(Vphy=Vphy),
                       family = bernoulli(), 
                       prior = prior2,
                       sample_prior = TRUE, chains = 4, cores =4, control=list(adapt_delta=0.9999), 
                       iter = 4500, warmup = 2500)

######

pp_check(mdnc_bern2,nsamples=100) 

np <- nuts_params(mdnc_bern2)
str(np)
# extract the number of divergence transitions
sum(subset(np, Parameter == "divergent__")$Value) #

summary(rhat(mdnc_bern2)) 

summary(mdnc_bern2,waic=FALSE, loo=FALSE)

plot(mdnc_bern2)

plot(conditional_effects(mdnc_bern2), points = FALSE) 

vcov(mdnc_bern2,correlations=TRUE)

save(mdnc_bern2,file="~/MacroEcology_group/members/Dylan/Hawaii_Invasion/Cleaned_Data/BRMS_PresAbs_MDWNC_sept2021.RData")

# identify non-significant interactions

load(file="~/MacroEcology_group/members/Dylan/Hawaii_Invasion/Cleaned_Data/BRMS_PresAbs_MDWNC_sept2021.RData")

get_variables(mdnc_bern2)[1:25]

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
# "b_AridInds:SM.wt.diff.lrrs"     "b_Soilz:SM.wt.diff.lrrs" "b_MATs:SM.wt.diff.lrrs"
# "b_AridInds:maxDBH.wt.diff.lrrs"  "MATs:maxDBH.wt.diff.lrrs"  "b_Soilz:maxDBH.wt.diff.lrrs" 
# "b_MATs:WD.wt.diff.lrrs" "b_HIIs:WD.wt.diff.lrrs" "b_Soilz:WD.wt.diff.lrrs" 
# "b_Soilz:PD.stands" ""b_MATs:PD.stands" ""b_AridInds:PD.stands"
# remove non-significant interactions

prior3 <- get_prior(PresAbs ~AridInds+MATs+HIIs+Soilz+ PD.stands+WD.wt.diff.lrrs+maxDBH.wt.diff.lrrs+SM.wt.diff.lrrs+
                      PD.stands:MATs+ PD.stands:AridInds+PD.stands:Soilz+
                      WD.wt.diff.lrrs:HIIs+MATs:WD.wt.diff.lrrs+Soilz:WD.wt.diff.lrrs+
                      AridInds:maxDBH.wt.diff.lrrs+MATs:maxDBH.wt.diff.lrrs+Soilz:maxDBH.wt.diff.lrrs+ 
                      SM.wt.diff.lrrs:MATs+AridInds:SM.wt.diff.lrrs+SM.wt.diff.lrrs:Soilz+ 
                      (1|Plot_Area2/PlotID)+ (1|Island)+(1|species)+ (1|gr(phylo, cov=Vphy)),
                    data = mdnc, 
                    data2=list(Vphy=Vphy),
                    family = bernoulli())

# fit model

mdnc_bern2_red <- brm(PresAbs ~ AridInds+MATs+HIIs+Soilz+ PD.stands+WD.wt.diff.lrrs+maxDBH.wt.diff.lrrs+SM.wt.diff.lrrs+
                        PD.stands:MATs+ PD.stands:AridInds+PD.stands:Soilz+
                        WD.wt.diff.lrrs:HIIs+MATs:WD.wt.diff.lrrs+Soilz:WD.wt.diff.lrrs+
                        AridInds:maxDBH.wt.diff.lrrs+MATs:maxDBH.wt.diff.lrrs+Soilz:maxDBH.wt.diff.lrrs+ 
                        SM.wt.diff.lrrs:MATs+AridInds:SM.wt.diff.lrrs+SM.wt.diff.lrrs:Soilz+
                        (1|Plot_Area2/PlotID)+ (1|Island)+(1|species)+ (1|gr(phylo, cov=Vphy)),
                      data = mdnc, 
                      data2=list(Vphy=Vphy),
                      family = bernoulli(),
                      prior=prior3,
                      sample_prior = TRUE, chains = 4, cores =4, control=list(adapt_delta=0.9999, max_treedepth=12), 
                      iter = 4500, warmup = 2500)

pp_check(mdnc_bern2_red,nsamples=100) 

np <- nuts_params(mdnc_bern2_red)
str(np)
# extract the number of divergence transitions
sum(subset(np, Parameter == "divergent__")$Value) #

summary(rhat(mdnc_bern2_red)) 

summary(mdnc_bern2_red,waic=FALSE, loo=FALSE)

save(mdnc_bern2,mdnc_bern2_red, file="~/MacroEcology_group/members/Dylan/Hawaii_Invasion/Cleaned_Data/BRMS_PresAbs_MDWNC_sept2021.RData")


# remove non-sig. interactions part II


prior4 <- get_prior(PresAbs ~AridInds+MATs+HIIs+Soilz+ PD.stands+WD.wt.diff.lrrs+maxDBH.wt.diff.lrrs+SM.wt.diff.lrrs+
                      PD.stands:MATs+ PD.stands:AridInds+
                      WD.wt.diff.lrrs:HIIs+Soilz:WD.wt.diff.lrrs+
                      AridInds:maxDBH.wt.diff.lrrs+MATs:maxDBH.wt.diff.lrrs+
                      SM.wt.diff.lrrs:MATs+AridInds:SM.wt.diff.lrrs+SM.wt.diff.lrrs:Soilz+ 
                      (1|Plot_Area2/PlotID)+ (1|Island)+(1|species)+ (1|gr(phylo, cov=Vphy)),
                    data = mdnc, 
                    data2=list(Vphy=Vphy),
                    family = bernoulli())

# fit model

mdnc_bern3_red <- brm(PresAbs ~ AridInds+MATs+HIIs+Soilz+ PD.stands+WD.wt.diff.lrrs+maxDBH.wt.diff.lrrs+SM.wt.diff.lrrs+
                        PD.stands:MATs+ PD.stands:AridInds+
                        WD.wt.diff.lrrs:HIIs+Soilz:WD.wt.diff.lrrs+
                        AridInds:maxDBH.wt.diff.lrrs+MATs:maxDBH.wt.diff.lrrs+
                        SM.wt.diff.lrrs:MATs+AridInds:SM.wt.diff.lrrs+SM.wt.diff.lrrs:Soilz+ 
                        (1|Plot_Area2/PlotID)+ (1|Island)+(1|species)+ (1|gr(phylo, cov=Vphy)),
                      data = mdnc, 
                      data2=list(Vphy=Vphy),
                      family = bernoulli(),
                      prior=prior4,
                      sample_prior = TRUE, chains = 4, cores =4, control=list(adapt_delta=0.9999, max_treedepth=12), 
                      iter = 4500, warmup = 2500)

pp_check(mdnc_bern3_red,nsamples=100) 

np <- nuts_params(mdnc_bern3_red)
str(np)
# extract the number of divergence transitions
sum(subset(np, Parameter == "divergent__")$Value) #

summary(rhat(mdnc_bern3_red)) 

summary(mdnc_bern3_red,waic=FALSE, loo=FALSE)

save(mdnc_bern2,mdnc_bern2_red, mdnc_bern3_red,file="~/MacroEcology_group/members/Dylan/Hawaii_Invasion/Cleaned_Data/BRMS_PresAbs_MDWNC_sept2021.RData")


# ##############
# # MDNN #######
# ##############
# 
# mdnn<-filter(local_all_paa,PDmeasure=="MDNN")
# mdnn$PDmeasure<-droplevels(mdnn$PDmeasure)
# #mdnn<-filter(mdnn, is.na(AGB)==FALSE)
# 
# #############
# # quick QC ##
# #############
# 
# length(unique(mdnn$PlotID)) # 460
# length(unique(mdnn$phylo)) #  43
# length(unique(tree$tip.label)) #  43
# 
# # transform + z-transform variables
# 
# mdnn$MAPs<-scale(mdnn$MAP,center=TRUE,scale=TRUE)
# mdnn$MATs<-scale(mdnn$MAT,center=TRUE, scale=TRUE)
# mdnn$Elev_ms<-scale(mdnn$Elev_m,center=TRUE, scale=TRUE)
# mdnn$AridInds<-scale(mdnn$AridInd,center=TRUE,scale=TRUE)
# mdnn$VPDs<-scale(mdnn$VPD,center=TRUE,scale=TRUE)
# mdnn$HIIs<-scale(mdnn$HII,center=TRUE,scale=TRUE)
# mdnn$HFPs<-scale(mdnn$HFP,center=TRUE,scale=TRUE)
# mdnn$SubstrateAge_range_lg<-log(mdnn$SubstrateAge_range)
# mdnn$Soilz<-scale(mdnn$SubstrateAge_range_lg,center=TRUE,scale=TRUE)
# mdnn$SeedWeights<-scale(log(mdnn$SeedWeight_g),center=TRUE,scale=TRUE)
# mdnn$WDs<-scale(mdnn$meanWD,center=TRUE,scale=TRUE)
# mdnn$maxHTs<-scale(mdnn$maxHT,center=TRUE,scale=TRUE)
# mdnn$PD.stands<-scale(mdnn$PD.stand,center=TRUE, scale=TRUE)
# 
# mdnn$S_PIE_lg<-log(mdnn$S_PIE)
# mdnn$S_PIEs<-scale(mdnn$S_PIE_lg,center=TRUE, scale=TRUE)
# mdnn$S_n_lg<-log(mdnn$S_n)
# mdnn$S_ns<-scale(mdnn$S_n_lg,center=TRUE, scale=TRUE)
# 
# mdnn$MPD_SES_lg<-log(mdnn$MPD_SES+0.00001)
# mdnn$MPD_SESs<-scale(mdnn$MPD_SES_lg,center=TRUE, scale=TRUE)
# mdnn$MNTD_SES_lg<-log(mdnn$MNTD_SES+0.00001)
# mdnn$MNTD_SESs<-scale(mdnn$MNTD_SES_lg,center=TRUE, scale=TRUE)
# 
# mdnn$eMPD_SES_lg<-log(mdnn$eMPD_SES+0.00001)
# mdnn$eMPD_SESs<-scale(mdnn$eMPD_SES_lg,center=TRUE, scale=TRUE)
# mdnn$eMNTD_SES_lg<-log(mdnn$eMNTD_SES+0.00001)
# mdnn$eMNTD_SESs<-scale(mdnn$eMNTD_SES_lg,center=TRUE, scale=TRUE)
# 
# mdnn$WD.diff.lrrs<-scale(mdnn$WD.diff.lrr,center=TRUE, scale=TRUE)
# mdnn$WD.wt.diff.lrrs<-scale(mdnn$WD.wt.diff.lrr,center=TRUE, scale=TRUE)
# mdnn$maxHT.diff.lrrs<-scale(mdnn$maxHT.diff.lrr,center=TRUE, scale=TRUE)
# mdnn$maxHT.wt.diff.lrrs<-scale(mdnn$maxHT.wt.diff.lrr,center=TRUE, scale=TRUE)
# mdnn$SM.diff.lrrs<-scale(mdnn$SM.diff.lrr,center=TRUE, scale=TRUE)
# mdnn$SM.wt.diff.lrrs<-scale(mdnn$SM.wt.diff.lrr,center=TRUE, scale=TRUE)
# 
# ##############################
# # corelations of predictors  #
# ##############################
# 
# mdnn<-ungroup(mdnn)
# mdnn_c<-select(mdnn, MAPs,MATs, VPDs,AridInds, HIIs, Soilz, PD.stands, Elev_ms, WDs, SeedWeights)
# 
# p_out<-corr.test(mdnn_c,method="pearson") # all below 0.7
# 
# # model structure:
# # PresAbs ~MAPs+PETs+ HIIs+Soilz+ AGBs+ PD.stands+WDs+maxDBHs+SeedWeights+
# #PD.stands:MAPs+ PD.stands:PETs+PD.stands:Soilz+PD.stands:HIIs+ PD.stands:AGBs+ AGBs:HIIs+
# 
# # get priors
# 
# # prior <- get_prior(PresAbs ~ MAPs+PETs+HIIs+Soilz+ PD.stands+WDs+SeedWeights+
# #                      PD.stands:MAPs+ PD.stands:PETs+PD.stands:Soilz+PD.stands:HIIs+
# #                      WDs:MAPs+WDs:PETs+WDs:Soilz+WDs:HIIs+
# #                      SeedWeights:MAPs+SeedWeights:PETs+SeedWeights:Soilz+SeedWeights:HIIs+ 
# #                      (1|Plot_Area2/PlotID)+ (1|phylo),data = mdnn, 
# #                    family = bernoulli())
# # 
# # # fit model
# # 
# # mdnn_bern <- brm(PresAbs ~MAPs+PETs+HIIs+Soilz+ PD.stands+WDs+SeedWeights+
# #                    PD.stands:MAPs+ PD.stands:PETs+PD.stands:Soilz+PD.stands:HIIs+
# #                    WDs:MAPs+WDs:PETs+WDs:Soilz+WDs:HIIs+
# #                    SeedWeights:MAPs+SeedWeights:PETs+SeedWeights:Soilz+SeedWeights:HIIs+ 
# #                    (1|Plot_Area2/PlotID)+(1|phylo), data = mdnn,
# #                  family = bernoulli(), cov_ranef = list(phylo = Vphy),
# #                  prior = prior,
# #                  sample_prior = TRUE, chains = 4, cores =4, 
# #                  iter = 4500, warmup = 2500) #,control=list(adapt_delta=0.95))
# 
# #mdnn_bernn <- kfold(mdnn_bern, K=10)
# 
# prior2 <- get_prior(PresAbs ~AridInds+MATs+HIIs+Soilz+ PD.stands+WD.wt.diff.lrrs+maxHT.wt.diff.lrrs+SM.wt.diff.lrrs+
#                       PD.stands:MATs+ PD.stands:AridInds+PD.stands:Soilz+PD.stands:HIIs+
#                       WD.wt.diff.lrrs:MATs+WD.wt.diff.lrrs:AridInds+WD.wt.diff.lrrs:Soilz+WD.wt.diff.lrrs:HIIs+
#                       maxHT.wt.diff.lrrs:MATs+maxHT.wt.diff.lrrs:AridInds+maxHT.wt.diff.lrrs:Soilz+maxHT.wt.diff.lrrs:HIIs+ 
#                       SM.wt.diff.lrrs:MATs+SM.wt.diff.lrrs:AridInds+SM.wt.diff.lrrs:Soilz+SM.wt.diff.lrrs:HIIs+ 
#                       (1|Plot_Area2/PlotID) +(1|phylo), data = mdnn,
#                     family = bernoulli())
# 
# mdnn_bern2 <- brm(PresAbs ~ AridInds+MATs+HIIs+Soilz+ PD.stands+WD.wt.diff.lrrs+maxHT.wt.diff.lrrs+SM.wt.diff.lrrs+
#                     PD.stands:MATs+ PD.stands:AridInds+PD.stands:Soilz+PD.stands:HIIs+
#                     WD.wt.diff.lrrs:MATs+WD.wt.diff.lrrs:AridInds+WD.wt.diff.lrrs:Soilz+WD.wt.diff.lrrs:HIIs+
#                     maxHT.wt.diff.lrrs:MATs+maxHT.wt.diff.lrrs:AridInds+maxHT.wt.diff.lrrs:Soilz+maxHT.wt.diff.lrrs:HIIs+ 
#                     SM.wt.diff.lrrs:MATs+SM.wt.diff.lrrs:AridInds+SM.wt.diff.lrrs:Soilz+SM.wt.diff.lrrs:HIIs+ 
#                     (1|Plot_Area2/PlotID)+ (1|phylo), data = mdnn,
#                  family = bernoulli(), cov_ranef = list(phylo = Vphy),
#                  prior = prior2,
#                  sample_prior = TRUE, chains = 4, cores =4,control=list(adapt_delta=0.99), 
#                  iter = 4500, warmup = 2500)
# 
# 
# #mdnn_bernn2 <- kfold(mdnn_bern2, K=10)
# 
# # compare random effect structures
# 
# # comp2<-loo_compare(mdnn_bernn,mdnn_bernn2) # A negative elpd_diff favors the first model.
# # print(comp2, digits=3)
# 
# ######
# 
# pp_check(mdnn_bern2,nsamples=1000)
# 
# np <- nuts_params(mdnn_bern2)
# str(np)
# # extract the number of divergence transitions
# sum(subset(np, Parameter == "divergent__")$Value) #
# 
# summary(rhat(mdnn_bern2)) 
# 
# summary(mdnn_bern2,waic=FALSE, loo=FALSE)
# 
# plot(mdnn_bern2)
# 
# #plot(marginal_effects(mdnc_lognormal), points = FALSE) 
# 
# vcov(mdnn_bern2,correlations=TRUE)
# 
# 
# save(mdnn_bern2, file="~/MacroEcology_group/members/Dylan/Hawaii_Invasion/Cleaned_Data/BRMS_Local_All_PresAbs_MDNN_output.RData")
# 
# #######################
# # remove interactions #
# #######################
# 
# 
# prior3 <- get_prior(PresAbs ~AridInds+MATs+HIIs+Soilz+ PD.stands+WD.wt.diff.lrrs+maxHT.wt.diff.lrrs+SM.wt.diff.lrrs+
#                       PD.stands:MATs+ PD.stands:AridInds+PD.stands:Soilz+
#                       WD.wt.diff.lrrs:HIIs+
#                       SM.wt.diff.lrrs:MATs+SM.wt.diff.lrrs:Soilz+ 
#                       (1|Plot_Area2/PlotID) +(1|phylo), data = mdnn,
#                     family = bernoulli())
# 
# mdnn_bern_red <- brm(PresAbs ~ AridInds+MATs+HIIs+Soilz+ PD.stands+WD.wt.diff.lrrs+maxHT.wt.diff.lrrs+SM.wt.diff.lrrs+
#                     PD.stands:MATs+ PD.stands:AridInds+PD.stands:Soilz+
#                     WD.wt.diff.lrrs:HIIs+
#                     SM.wt.diff.lrrs:MATs+SM.wt.diff.lrrs:Soilz+ 
#                     (1|Plot_Area2/PlotID)+ (1|phylo), data = mdnn,
#                   family = bernoulli(), cov_ranef = list(phylo = Vphy),
#                   prior = prior3,
#                   sample_prior = TRUE, chains = 4, cores =4,control=list(adapt_delta=0.99), 
#                   iter = 4500, warmup = 2500)
# 
# ######
# 
# pp_check(mdnn_bern_red,nsamples=1000)
# 
# np <- nuts_params(mdnn_bern_red)
# str(np)
# # extract the number of divergence transitions
# sum(subset(np, Parameter == "divergent__")$Value) #
# 
# summary(rhat(mdnn_bern_red)) 
# 
# summary(mdnn_bern_red,waic=FALSE, loo=FALSE)
# 
# plot(mdnn_bern_red)
# 
# #plot(marginal_effects(mdnc_lognormal), points = FALSE) 
# 
# vcov(mdnn_bern_red,correlations=TRUE)
# 
# save(mdnn_bern2,mdnn_bern_red,
#      file="~/MacroEcology_group/members/Dylan/Hawaii_Invasion/Cleaned_Data/BRMS_Local_All_PresAbs_MDNN_output.RData")