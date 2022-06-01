################################
# Hierarchical Bayesian Model ##
# Response Variable: Rel.Abund.#
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

load("Cleaned_Data/Local_All_RelAbund_PresAbs_Final.RData")

local_all_raa$PlotID<-as.factor(local_all_raa$PlotID)

local_all_raa$Plot_Area2<-as.factor(local_all_raa$Plot_Area)

colnames(local_all_raa)[17]<-"phylo"

local_all_raa$phylo<-as.character(local_all_raa$phylo)

local_all_raa$species<-local_all_raa$phylo

######################
# add in phylogeny   #
######################

# tree<-read.tree("Cleaned_Data/Local_Aliens_SizeClass_All_phylogeny.tre")

tree<-local_tree_raa

#scale co-variance matrix

nspp<-length(unique(tree$tip.label))
Vphy <- vcv(tree)
Vphy <- Vphy/(det(Vphy)^(1/nspp))

###############
# MDNC ########
###############

mdnc<-filter(local_all_raa,PDmeasure=="MDWNC")
mdnc$PDmeasure<-droplevels(mdnc$PDmeasure)
# mdnc<-filter(mdnc, is.na(AGB)==FALSE)

#############
# quick QC ##
#############

length(unique(mdnc$PlotID)) # 167
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
mdnc$RelAbund<-mdnc$Abundance_ha/mdnc$Total_Abundance_ha

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
mdnc_c<-dplyr::select(mdnc, MATs, AridInds,HIIs, Soilz,
                PD.stands, WD.wt.diff.lrrs,maxDBH.wt.diff.lrrs,
                SM.wt.diff.lrrs)

p_out<-corr.test(mdnc_c,method="pearson") # all below 0.7

# below 0.7

test <- lmer(RelAbund ~  AridInds+MATs+HIIs+Soilz+ PD.stands+WD.wt.diff.lrrs+maxDBH.wt.diff.lrrs+SM.wt.diff.lrrs+
                     PD.stands:MATs+ PD.stands:AridInds+PD.stands:Soilz+PD.stands:HIIs+
                     WD.wt.diff.lrrs:MATs+WD.wt.diff.lrrs:AridInds+WD.wt.diff.lrrs:Soilz+WD.wt.diff.lrrs:HIIs+
                     maxDBH.wt.diff.lrrs:MATs+maxDBH.wt.diff.lrrs:AridInds+maxDBH.wt.diff.lrrs:Soilz+maxDBH.wt.diff.lrrs:HIIs+ 
                     SM.wt.diff.lrrs:MATs+SM.wt.diff.lrrs:AridInds+SM.wt.diff.lrrs:Soilz+SM.wt.diff.lrrs:HIIs+ 
                     (1|PlotID)+ (1|phylo),data = mdnc) 

vif_mod<-vif(test) 

#########################
# STEP ONE ##############
# 1. select RE          #
# FULL MODEL            #
#########################
## log normal or beta
# get priors

prior <- get_prior(RelAbund ~  AridInds+MATs+HIIs+Soilz+ PD.stands+WD.wt.diff.lrrs+maxDBH.wt.diff.lrrs+SM.wt.diff.lrrs+
                     PD.stands:MATs+ PD.stands:AridInds+PD.stands:Soilz+PD.stands:HIIs+
                     WD.wt.diff.lrrs:MATs+WD.wt.diff.lrrs:AridInds+WD.wt.diff.lrrs:Soilz+WD.wt.diff.lrrs:HIIs+
                     maxDBH.wt.diff.lrrs:MATs+maxDBH.wt.diff.lrrs:AridInds+maxDBH.wt.diff.lrrs:Soilz+maxDBH.wt.diff.lrrs:HIIs+ 
                     SM.wt.diff.lrrs:MATs+SM.wt.diff.lrrs:AridInds+SM.wt.diff.lrrs:Soilz+SM.wt.diff.lrrs:HIIs+ 
                     (1|Plot_Area2/PlotID)+  (1|Island)+ (1|species)+ (1|gr(phylo, cov=Vphy)),
                   data = mdnc, 
                   data2=list(Vphy=Vphy),
                   family = Beta())

# fit beta model

mdnc_beta <- brm(RelAbund ~ AridInds+MATs+HIIs+Soilz+ PD.stands+WD.wt.diff.lrrs+maxDBH.wt.diff.lrrs+SM.wt.diff.lrrs+
                        PD.stands:MATs+ PD.stands:AridInds+PD.stands:Soilz+PD.stands:HIIs+
                        WD.wt.diff.lrrs:MATs+WD.wt.diff.lrrs:AridInds+WD.wt.diff.lrrs:Soilz+WD.wt.diff.lrrs:HIIs+
                        maxDBH.wt.diff.lrrs:MATs+maxDBH.wt.diff.lrrs:AridInds+maxDBH.wt.diff.lrrs:Soilz+maxDBH.wt.diff.lrrs:HIIs+ 
                        SM.wt.diff.lrrs:MATs+SM.wt.diff.lrrs:AridInds+SM.wt.diff.lrrs:Soilz+SM.wt.diff.lrrs:HIIs+ 
                        (1|Plot_Area2/PlotID)+ (1|Island)+(1|species)+ (1|gr(phylo, cov=Vphy)),
                      data = mdnc, 
                      data2=list(Vphy=Vphy), family = Beta(),
                 prior = prior,
                 sample_prior = TRUE, chains = 4, cores =4,  
                 control=list(adapt_delta=0.99), 
                 iter = 4500, warmup = 2500)

mdnc_beta_kfold <- kfold(mdnc_beta, chains = 4, cores = 4)

#  fit lognormal model

prior2 <- get_prior(RelAbund ~  AridInds+MATs+HIIs+Soilz+ PD.stands+WD.wt.diff.lrrs+maxDBH.wt.diff.lrrs+SM.wt.diff.lrrs+
                     PD.stands:MATs+ PD.stands:AridInds+PD.stands:Soilz+PD.stands:HIIs+
                     WD.wt.diff.lrrs:MATs+WD.wt.diff.lrrs:AridInds+WD.wt.diff.lrrs:Soilz+WD.wt.diff.lrrs:HIIs+
                     maxDBH.wt.diff.lrrs:MATs+maxDBH.wt.diff.lrrs:AridInds+maxDBH.wt.diff.lrrs:Soilz+maxDBH.wt.diff.lrrs:HIIs+ 
                     SM.wt.diff.lrrs:MATs+SM.wt.diff.lrrs:AridInds+SM.wt.diff.lrrs:Soilz+SM.wt.diff.lrrs:HIIs+ 
                     (1|Plot_Area2/PlotID)+ (1|Island)+(1|species)+ (1|gr(phylo, cov=Vphy)),
                   data = mdnc, 
                   data2=list(Vphy=Vphy),
                   family = lognormal())

# fit model

mdnc_lognormal <- brm(RelAbund ~ AridInds+MATs+HIIs+Soilz+ PD.stands+WD.wt.diff.lrrs+maxDBH.wt.diff.lrrs+SM.wt.diff.lrrs+
                   PD.stands:MATs+ PD.stands:AridInds+PD.stands:Soilz+PD.stands:HIIs+
                   WD.wt.diff.lrrs:MATs+WD.wt.diff.lrrs:AridInds+WD.wt.diff.lrrs:Soilz+WD.wt.diff.lrrs:HIIs+
                   maxDBH.wt.diff.lrrs:MATs+maxDBH.wt.diff.lrrs:AridInds+maxDBH.wt.diff.lrrs:Soilz+maxDBH.wt.diff.lrrs:HIIs+ 
                   SM.wt.diff.lrrs:MATs+SM.wt.diff.lrrs:AridInds+SM.wt.diff.lrrs:Soilz+SM.wt.diff.lrrs:HIIs+ 
                   (1|Plot_Area2/PlotID)+ (1|Island)+(1|species)+ (1|gr(phylo, cov=Vphy)),
                 data = mdnc, 
                 data2=list(Vphy=Vphy), family = lognormal(),
                 prior = prior2,
                 sample_prior = TRUE, chains = 4, cores =4,  
                 control=list(adapt_delta=0.999), 
                 iter = 4500, warmup = 2500)

mdnc_lognormal_kfold <- kfold(mdnc_lognormal, chains = 4, cores = 4)

# compare
comp<-loo_compare(mdnc_beta_kfold, mdnc_lognormal_kfold) # A negative elpd_diff favors the first model.


# keep mdnc_lognormal2 b/c it accounts for plot area

############################ 
# Step 2                 ##
# 3. summarize         #####
############################

pp_check( mdnc_lognormal,nsamples=100) +scale_x_continuous(trans="log")

np <- nuts_params(mdnc_lognormal)
str(np)
# extract the number of divergence transitions
sum(subset(np, Parameter == "divergent__")$Value) #

summary(rhat(mdnc_lognormal)) 

summary(mdnc_lognormal,waic=FALSE, loo=FALSE)

plot(mdnc_lognormal)

plot(conditional_effects(mdnc_lognormal), points = FALSE) 

vcov(mdnc_lognormal,correlations=TRUE)

save(mdnc_lognormal,comp, file="~/MacroEcology_group/members/Dylan/Hawaii_Invasion/Cleaned_Data/BRMS_RelAbund_MDWNC_sept2021.RData")

#-------------------------------------#
#  test [remove non-sig interactions] #
#------------------------------------ #

# load(file="~/MacroEcology_group/members/Dylan/Hawaii_Invasion/Cleaned_Data/BRMS_Local_All_RelAbund_MDWNC_output.RData")

get_variables(mdnc_lognormal)[1:25]

mdnc_coef2<-mdnc_lognormal %>%
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

mdnc_coef2<-dplyr::filter(mdnc_coef2, terms!="b_Intercept")

mdnc_coef2_pp<-ggplot(mdnc_coef2, aes(y=terms, x=value))+
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

#no interactions

#  plot nested in area

prior3 <- get_prior(RelAbund ~  AridInds+MATs+HIIs+Soilz+ PD.stands+WD.wt.diff.lrrs+maxDBH.wt.diff.lrrs+SM.wt.diff.lrrs+
                      (1|Plot_Area2/PlotID)+ (1|Island)+(1|species)+ (1|gr(phylo, cov=Vphy)),
                    data = mdnc, 
                    data2=list(Vphy=Vphy),
                    family = lognormal())

# fit model

mdnc_lognormal_red <- brm(RelAbund ~  AridInds+MATs+HIIs+Soilz+ PD.stands+WD.wt.diff.lrrs+maxDBH.wt.diff.lrrs+
                            SM.wt.diff.lrrs+
                            (1|Plot_Area2/PlotID)+ (1|Island)+(1|species)+ (1|gr(phylo, cov=Vphy)),
                          data = mdnc, 
                          data2=list(Vphy=Vphy), family = lognormal(),
                          prior = prior3,
                          sample_prior = TRUE, chains = 4, cores =4,  
                          control=list(adapt_delta=0.99), 
                          iter = 4500, warmup = 2500)

pp_check( mdnc_lognormal_red,nsamples=100) +scale_x_continuous(trans="log")

np <- nuts_params(mdnc_lognormal_red)
str(np)
# extract the number of divergence transitions
sum(subset(np, Parameter == "divergent__")$Value) #

summary(rhat(mdnc_lognormal_red)) 

summary(mdnc_lognormal_red,waic=FALSE, loo=FALSE)

plot(mdnc_lognormal_red)

plot(conditional_effects(mdnc_lognormal_red), points = FALSE) 

save(mdnc_lognormal,comp,mdnc_lognormal_red, file="~/MacroEcology_group/members/Dylan/Hawaii_Invasion/Cleaned_Data/BRMS_RelAbund_MDWNC_sept2021.RData")

# ##############
# # MDNN #######
# ##############
# 
# mdnn<-filter(local_all_raa,PDmeasure=="MDNN")
# 
# mdnn$PDmeasure<-droplevels(mdnn$PDmeasure)
# # mdnn<-filter(mdnn, is.na(AGB)==FALSE)
# 
# #############
# # quick QC ##
# #############
# 
# length(unique(mdnn$PlotID)) # 170
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
# mdnn<-ungroup(mdnn)
# mdnn_c<-select(mdnn, MAPs,MATs, AridInds, VPDs,HIIs, Soilz, PD.stands, Elev_ms, WDs, SeedWeights)
# 
# p_out<-corr.test(mdnn_c,method="pearson") # all below 0.7
# 
# # model structure:
# # PresAbs ~MAPs+PETs+ HIIs+Soilz+ AGBs+ PD.stands+WDs+maxDBHs+SeedWeights+
# #PD.stands:MAPs+ PD.stands:PETs+PD.stands:Soilz+PD.stands:HIIs+ PD.stands:AGBs+ AGBs:HIIs+
# 
# # get priors
# 
# #  plot nested in area
# 
# prior2 <- get_prior(RelAbund ~ AridInds+MATs+HIIs+Soilz+ PD.stands+WD.wt.diff.lrrs+maxHT.wt.diff.lrrs+SM.wt.diff.lrrs+
#                       PD.stands:MATs+ PD.stands:AridInds+PD.stands:Soilz+PD.stands:HIIs+
#                       WD.wt.diff.lrrs:MATs+WD.wt.diff.lrrs:AridInds+WD.wt.diff.lrrs:Soilz+WD.wt.diff.lrrs:HIIs+
#                       maxHT.wt.diff.lrrs:MATs+maxHT.wt.diff.lrrs:AridInds+maxHT.wt.diff.lrrs:Soilz+maxHT.wt.diff.lrrs:HIIs+ 
#                       SM.wt.diff.lrrs:MATs+SM.wt.diff.lrrs:AridInds+SM.wt.diff.lrrs:Soilz+SM.wt.diff.lrrs:HIIs+ 
#                       (1|Plot_Area2/PlotID)+ (1|phylo),data = mdnn, 
#                     family = lognormal())
# 
# # fit model
# 
# mdnn_lognormal2 <- brm(RelAbund ~ AridInds+MATs+HIIs+Soilz+ PD.stands+WD.wt.diff.lrrs+maxHT.wt.diff.lrrs+SM.wt.diff.lrrs+
#                          PD.stands:MATs+ PD.stands:AridInds+PD.stands:Soilz+PD.stands:HIIs+
#                          WD.wt.diff.lrrs:MATs+WD.wt.diff.lrrs:AridInds+WD.wt.diff.lrrs:Soilz+WD.wt.diff.lrrs:HIIs+
#                          maxHT.wt.diff.lrrs:MATs+maxHT.wt.diff.lrrs:AridInds+maxHT.wt.diff.lrrs:Soilz+maxHT.wt.diff.lrrs:HIIs+ 
#                          SM.wt.diff.lrrs:MATs+SM.wt.diff.lrrs:AridInds+SM.wt.diff.lrrs:Soilz+SM.wt.diff.lrrs:HIIs+ 
#                          (1|Plot_Area2/PlotID)+ (1|phylo), data = mdnn,
#                        family = lognormal(), cov_ranef = list(phylo = Vphy),
#                        prior = prior2,
#                        sample_prior = TRUE, chains = 4, cores =4, control=list(adapt_delta=0.99), 
#                        iter = 4500, warmup = 2500)
# 
# 
# ######
# 
# pp_check(mdnn_lognormal2,nsamples=1000) +scale_x_continuous(trans="log")
# 
# np <- nuts_params(mdnn_lognormal2)
# str(np)
# # extract the number of divergence transitions
# sum(subset(np, Parameter == "divergent__")$Value) #
# 
# summary(rhat(mdnn_lognormal2)) 
# 
# summary(mdnn_lognormal2,waic=FALSE, loo=FALSE)
# 
# plot(mdnn_lognormal2)
# 
# plot(marginal_effects(mdnn_lognormal2), points = FALSE) 
# 
# vcov(mdnn_lognormal2,correlations=TRUE)
# 
# # re-fit with only significant interactions
# 
# 
# prior3 <- get_prior(RelAbund ~ AridInds+MATs+HIIs+Soilz+ PD.stands+WD.wt.diff.lrrs+maxHT.wt.diff.lrrs+SM.wt.diff.lrrs+
#                       WD.wt.diff.lrrs:MATs+
#                       maxHT.wt.diff.lrrs:AridInds+ 
#                       SM.wt.diff.lrrs:HIIs+ 
#                       (1|Plot_Area2/PlotID)+ (1|phylo),data = mdnn, 
#                     family = lognormal())
# 
# # fit model
# 
# mdnn_lognormal_red <- brm(RelAbund ~ AridInds+MATs+HIIs+Soilz+ PD.stands+WD.wt.diff.lrrs+maxHT.wt.diff.lrrs+SM.wt.diff.lrrs+
#                          WD.wt.diff.lrrs:MATs+
#                          maxHT.wt.diff.lrrs:AridInds+ 
#                          SM.wt.diff.lrrs:HIIs+ 
#                          (1|Plot_Area2/PlotID)+ (1|phylo), data = mdnn,
#                        family = lognormal(), cov_ranef = list(phylo = Vphy),
#                        prior = prior3,
#                        sample_prior = TRUE, chains = 4, cores =4, control=list(adapt_delta=0.99), 
#                        iter = 4500, warmup = 2500)
# 
# ####
# save(mdnn_lognormal2,mdnn_lognormal_red, file="~/MacroEcology_group/members/Dylan/Hawaii_Invasion/Cleaned_Data/BRMS_Local_All_RelAbund_MDNN_output.RData")