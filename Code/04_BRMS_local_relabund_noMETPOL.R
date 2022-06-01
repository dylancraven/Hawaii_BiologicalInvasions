################################
# Hierarchical Bayesian Model ##
# Response Variable: Rel.Abund.#
# Scale: local                 #
# Size Class: all              #
################################
# excluding METPOL #############
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

load("~/MacroEcology_group/members/Dylan/Hawaii_Invasion/Cleaned_Data/Local_All_RelAbund_PresAbs_NoMETPOL_Final.RData")

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
#mdnc$PDmeasure<-droplevels(mdnc$PDmeasure)
# mdnc<-filter(mdnc, is.na(AGB)==FALSE)

#############
# quick QC ##
#############

length(unique(mdnc$PlotID)) # 112
length(unique(mdnc$phylo)) #  36
length(unique(tree$tip.label)) #  36

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

#  plot nested in area

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
                      control=list(adapt_delta=0.99), 
                      iter = 4500, warmup = 2500)

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

save(mdnc_lognormal,file="Cleaned_Data/BRMS_RelAbund_MDWNC_noMETPOL_sept2021.RData")

#-------------------------------------#
#  test [remove non-sig interactions] #
#------------------------------------ #

load(file="~/MacroEcology_group/members/Dylan/Hawaii_Invasion/Cleaned_Data/BRMS_RelAbund_MDWNC_noMETPOL_sept2021.RData")

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


# second removed all interactions

#  plot nested in area

prior3 <- get_prior(RelAbund ~ AridInds+MATs+HIIs+Soilz+ PD.stands+WD.wt.diff.lrrs+maxDBH.wt.diff.lrrs+SM.wt.diff.lrrs+
                      (1|Plot_Area2/PlotID)+ (1|Island)+(1|species)+ (1|gr(phylo, cov=Vphy)),
                    data = mdnc, 
                    data2=list(Vphy=Vphy),
                    family = lognormal())

# fit model

mdnc_lognormal2_red <- brm(RelAbund ~  AridInds+MATs+HIIs+Soilz+ PD.stands+WD.wt.diff.lrrs+maxDBH.wt.diff.lrrs+SM.wt.diff.lrrs+
                             (1|Plot_Area2/PlotID)+ (1|Island)+(1|species)+ (1|gr(phylo, cov=Vphy)),
                           data = mdnc, 
                           data2=list(Vphy=Vphy), family = lognormal(),
                           prior = prior3,
                           sample_prior = TRUE, chains = 4, cores =4,  
                           control=list(adapt_delta=0.99), 
                           iter = 4500, warmup = 2500)

pp_check( mdnc_lognormal2_red,nsamples=100) +scale_x_continuous(trans="log")

np <- nuts_params(mdnc_lognormal2_red)
str(np)
# extract the number of divergence transitions
sum(subset(np, Parameter == "divergent__")$Value) #

summary(rhat(mdnc_lognormal2_red)) 

summary(mdnc_lognormal2_red,waic=FALSE, loo=FALSE)

save(mdnc_lognormal,mdnc_lognormal2_red, file="~/MacroEcology_group/members/Dylan/Hawaii_Invasion/Cleaned_Data/BRMS_RelAbund_MDWNC_noMETPOL_sept2021.RData")

