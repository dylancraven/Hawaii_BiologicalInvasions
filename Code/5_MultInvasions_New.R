#----------------#
# Alien SppN or  #
# Relative Inv.  #
#----------------#

require(tidyverse)
require(brms)
require(sp)
require(spdep) 

#--data fun------# 

plot<-read.csv("Cleaned_Data/HawIslands_RelAbund_PresAbs_Plot.csv",sep=",",header=T)

local_sel<-plot %>% 
        dplyr::filter(., Native_Status=="native") %>% 
        group_by(PlotID) %>% 
        summarise(Abund_Tot=sum(Abundance_ha)) %>% 
        dplyr::filter(.,Abund_Tot>0)

local_sel<-unique(local_sel$PlotID)

plot<-dplyr::filter(plot, PlotID %in% local_sel)

env_1km<-read.csv("Cleaned_Data/Hawaii_newclimate_1km.csv",sep=",",header = T)

env_1km<-dplyr::select(env_1km, PlotID=PlotIDn,Long_Dec, Lat_Dec, 
                       MAT=newMAT,MAP=newMAP, Elev_m= newElev, VPD=newVPD, PET=newPET, 
                       AridInd=newAridInd, HFP, HII,SubstrateAge_range,SubstrateAge_code)%>%
                       distinct(.)

plot_size<-plot %>% 
           select(., PlotID, Plot_Area) %>% 
           distinct(.)

#--Abundance / Spp Richness------# 

alien_native_abund<-plot%>%
        group_by(Island,PlotID,Native_Status)%>%
        summarize(Abundance=sum(Abundance_ha))%>%
        pivot_wider(id_cols=c(Island,PlotID), names_from =c(Native_Status),
                    values_from = c(Abundance))%>%
        mutate( Abundance_both=alien+native,
                AlienAbund_Prop=alien/Abundance_both)%>%
        left_join(.,plot_size, by="PlotID") %>% 
        select(.,Island, PlotID,Plot_Area, Alien_abundance=alien, Native_abundance=native,Abundance_both,AlienAbund_Prop)

alien_native_sppn<-plot%>%
        filter(Abundance_ha>0)%>%
        group_by(Island,PlotID,Native_Status)%>%
        summarize(SppN=length(unique(SPP_CODE3A)))%>%
        pivot_wider(id_cols=c(Island,PlotID), names_from =c(Native_Status),
                    values_from = c(SppN))%>%
        mutate( alien=ifelse(is.na(alien)==TRUE,0,alien),
                SppN_both=alien+native,
                AlienSppN_Prop=alien/SppN_both)%>%
        select(.,Island, PlotID,Alien_SppN=alien, Native_SppN=native,SppN_both,AlienSppN_Prop)

alien_native<-alien_native_abund%>%
        left_join(.,alien_native_sppn, by=c("Island","PlotID"))%>%
        left_join(., env_1km, by="PlotID")%>%
        drop_na(.)

save(alien_native, file="Cleaned_Data/MultiInv_dat.RData")

#--------------------------#
#------ make models -------#
#--------------------------#

load(file="Cleaned_Data/MultiInv_dat.RData")

# scale predictors
alien_native$MAPs<-scale(alien_native$MAP,center=TRUE,scale=TRUE)
alien_native$MATs<-scale(alien_native$MAT,center=TRUE, scale=TRUE)
alien_native$Elev_ms<-scale(alien_native$Elev_m,center=TRUE, scale=TRUE)
alien_native$AridInd<-alien_native$AridInd*-1
alien_native$AridInds<-scale(alien_native$AridInd,center=TRUE,scale=TRUE)
alien_native$VPDs<-scale(alien_native$VPD,center=TRUE,scale=TRUE)
alien_native$HIIs<-scale(alien_native$HII,center=TRUE,scale=TRUE)
alien_native$HFPs<-scale(alien_native$HFP,center=TRUE,scale=TRUE)
alien_native$SubstrateAge_range_lg<-log(alien_native$SubstrateAge_range)
alien_native$Soilz<-scale(alien_native$SubstrateAge_range_lg,center=TRUE,scale=TRUE)
alien_native$Plot_Area<-as.factor(alien_native$Plot_Area)

alien_native<-ungroup(alien_native)

alien_native<-filter(alien_native, Alien_abundance>0)

#-----------#
# SppN Prop #
#-----------#

prior<- get_prior(AlienSppN_Prop~ AridInds+MATs+HIIs+Soilz+
                             (1|Island)+(1|Plot_Area),data = alien_native, 
                     family = Beta())

SppN_prop <- brm(AlienSppN_Prop~ AridInds+MATs+HIIs+Soilz+
                           (1|Island)+(1|Plot_Area),data = alien_native, 
                   prior = prior,
                   family = Beta(),
                   sample_prior = TRUE, chains = 4, cores =4,
                   iter = 4500, warmup = 2500,
                   control=list(adapt_delta=0.99999))

# kfold1 <- kfold(SppN_prop, chains = 4, cores = 4)
# 
# prior_b<- get_prior(AlienSppN_Prop~ AridInds+MATs+HIIs+Soilz+
#                           (1|Island),data = alien_native, 
#                   family = gaussian())
# 
# SppN_prop_b <- brm(AlienSppN_Prop~ AridInds+MATs+HIIs+Soilz+
#                          (1|Island),data = alien_native, 
#                  prior = prior_b,
#                  family = gaussian(),
#                  sample_prior = TRUE, chains = 4, cores =4,
#                  iter = 4500, warmup = 2500,
#                  control=list(adapt_delta=0.9999))
# 
# kfold2 <- kfold(SppN_prop_b, chains = 4, cores = 4)
# 
# gaus_v_beta<-loo_compare(kfold1, kfold2) # 

#--model check--# 

pp_check(SppN_prop,nsamples=100)

np <- nuts_params(SppN_prop)
str(np)
# extract the number of divergence transitions
sum(subset(np, Parameter == "divergent__")$Value) #

summary(rhat(SppN_prop)) 

plot(SppN_prop)

summary(SppN_prop,waic=TRUE, loo=TRUE)

bayes_R2(SppN_prop)

save(SppN_prop, 
     file="Cleaned_Data/Haw_MultInvasions.RData")

# make  spatial term (not necessary)
# load(file="Cleaned_Data/Haw_MultInvasions.RData")
# coords <- cbind(alien_native$Long_Dec, alien_native$Lat_Dec)
# myDist = 25
# res<-data.frame(residuals(SppN_prop, summary = TRUE))
#  res<-res$Estimate
#  rac <- autocov_dist(res, coords, nbs = myDist,
#                      type = "inverse", zero.policy = TRUE, style = "W", longlat=T) # longlat=T interprets myDist = 300 as 300 km.
# 
# alien_native$rac<-rac
# 
# # re-fit model
# 
# prior_s <- get_prior(Alien_SppN|trials(SppN_both)~ AridInds+MATs+HIIs+Soilz+rac+
#                              (1|Island),data = alien_native, 
#                      family = binomial())
# 
# SppN_prop_s <- brm(Alien_SppN|trials(SppN_both)~ AridInds+MATs+HIIs+Soilz+rac+
#                            (1|Island),data = alien_native, 
#                  prior = prior_s,
#                  family = binomial(),
#                  sample_prior = TRUE, chains = 4, cores =4,
#                  iter = 4500, warmup = 2500,
#                  control=list(adapt_delta=0.99))
# 
# # #--model check--# 
#  
# pp_check(SppN_prop_s,nsamples=100)
#  
# np <- nuts_params(SppN_prop_s)
# str(np)
# # # extract the number of divergence transitions
# sum(subset(np, Parameter == "divergent__")$Value) #
# 
# summary(rhat(SppN_prop_s))
# # 
# summary(SppN_prop_s,waic=TRUE, loo=TRUE)
# # 
# plot(SppN_prop_s)
# # 
# bayes_R2(SppN_prop_s)
# # 
# save(SppN_prop,SppN_prop_s,
#      file="Cleaned_Data/Haw_MultInvasions.RData")

#-----------------#
# AlienAbund_Prop #
#-----------------#

prior2 <- get_prior(AlienAbund_Prop  ~ AridInds+MATs+HIIs+Soilz+
                     (1|Island)+(1|Plot_Area),data = alien_native, 
                    family = Beta())

Abund_prop <- brm(AlienAbund_Prop  ~ AridInds+MATs+HIIs+Soilz+
                   (1|Island)+(1|Plot_Area),data = alien_native, 
                 prior = prior2,
                 family = Beta(),
                 sample_prior = TRUE, chains = 4, cores =4,
                 iter = 4500, warmup = 2500,
                 control=list(adapt_delta=0.99))

#--model check--# 

pp_check(Abund_prop,nsamples=100)

np <- nuts_params(Abund_prop)
str(np)
# extract the number of divergence transitions
sum(subset(np, Parameter == "divergent__")$Value) #

summary(rhat(Abund_prop))  

summary(Abund_prop,waic=TRUE, loo=TRUE)

plot(Abund_prop)

bayes_R2(Abund_prop)

#load(file="Cleaned_Data/Haw_MultInvasions.RData")
save(SppN_prop, Abund_prop,
     file="Cleaned_Data/Haw_MultInvasions.RData")

# make  spatial term (not needed)

load(file="Cleaned_Data/Haw_MultInvasions.RData")
coords <- cbind(alien_native$Long_Dec, alien_native$Lat_Dec)
myDist = 25
res2<-data.frame(residuals(Abund_prop, summary = TRUE))
res2<-res2$Estimate
rac2 <- autocov_dist(res2, coords, nbs = myDist,
                    type = "inverse", zero.policy = TRUE, style = "W", longlat=T) # longlat=T interprets myDist = 300 as 300 km.

alien_native$rac2<-rac2

prior3 <- get_prior(AlienAbund_Prop ~ AridInds+MATs+HIIs+Soilz+rac2+
                      (1|Island),data = alien_native, 
                    family = zero_one_inflated_beta())

Abund_prop2 <- brm(AlienAbund_Prop ~ AridInds+MATs+HIIs+Soilz+rac2+
                    (1|Island),data = alien_native, 
                  prior = prior3,
                  family = zero_one_inflated_beta(),
                  sample_prior = TRUE, chains = 4, cores =4,
                  iter = 4500, warmup = 2500,
                  control=list(adapt_delta=0.9999,max_treedepth=12))

#--model check--# 

pp_check(Abund_prop2,nsamples=100)

np <- nuts_params(Abund_prop2)
str(np)
# extract the number of divergence transitions
sum(subset(np, Parameter == "divergent__")$Value) #

summary(rhat(Abund_prop2)) 

summary(Abund_prop2,waic=TRUE, loo=TRUE)

plot(Abund_prop2)

bayes_R2(Abund_prop2)

save(SppN_prop,SppN_prop_s,Abund_prop,Abund_prop2,
     file="Cleaned_Data/Haw_MultInvasions.RData")
