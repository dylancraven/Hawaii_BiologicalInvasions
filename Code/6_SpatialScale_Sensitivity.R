#----------------#
# Spatial bias?  #
# species-level  #
# occurrence and #
# abundance      #
#----------------#

require(tidyverse)
require(brms)
require(cowplot)

load("Cleaned_Data/BRMS_PresAbs_MDWNC_sept2021.RData")
rm("mdnc_bern2","mdnc_bern2_red")

load("Cleaned_Data/BRMS_RelAbund_MDWNC_sept2021.RData")
rm(mdnc_lognormal)

plot_area_n<-mdnc_bern3_red$data %>% 
  mutate(Plot_Area2=as.character(Plot_Area2)) %>% 
  dplyr::group_by(Plot_Area2) %>% 
  dplyr::summarize(N=length(Plot_Area2)) %>% 
  mutate(Plot_Area2=round(as.numeric(Plot_Area2), digits=2))

#-------------------#
# Estimate accuracy #
# pres/abs.         #
#-------------------#

PA_accuracy <- predict(mdnc_bern3_red) %>% 
  data.frame(.) %>% 
  bind_cols(., mdnc_bern3_red$data) %>% 
  mutate(Accuracy=(Estimate-0.5)/0.5,
         Plot_Area2=as.character(Plot_Area2)) %>%
  group_by(Plot_Area2) %>% 
  summarise(ci = list(mean_cl_normal(Accuracy) %>% 
                        rename(accurancy_xbar=y, accuracy_lwr=ymin, accuracy_upr=ymax))) %>% 
  unnest(cols = c(ci)) %>% 
  mutate(Plot_Area2=round(as.numeric(Plot_Area2),digits=2))

#-------------------#
# Estimate precision#
# pres/abs.         #
#-------------------#

PA_precision <- data.frame(ranef(mdnc_bern3_red)$Plot_Area2) %>% 
  rownames_to_column(.,var="Plot_Area2") %>%  
  mutate(Precision=1/Est.Error.Intercept,
         Plot_Area2=as.character(Plot_Area2)) %>%
  mutate(Plot_Area2=round(as.numeric(Plot_Area2),digits=2)) %>% 
  select(., Plot_Area2, Estimate.Intercept,Q2.5.Intercept,Q97.5.Intercept,Precision) %>% 
  left_join(.,plot_area_n,by="Plot_Area2")

#-------------------#
# Estimate accuracy #
# dominance         #
#-------------------#

RelAbund_accuracy <- predict(mdnc_lognormal_red) %>% 
  data.frame(.) %>% 
  bind_cols(., mdnc_lognormal_red$data) %>% 
  mutate(Accuracy=(Estimate-RelAbund)/RelAbund,
         Plot_Area2=as.character(Plot_Area2)) %>%
  group_by(Plot_Area2) %>% 
  summarise(ci = list(mean_cl_normal(Accuracy) %>% 
                        rename(accurancy_xbar=y, accuracy_lwr=ymin, accuracy_upr=ymax))) %>% 
  unnest(cols = c(ci)) %>% 
  mutate(Plot_Area2=round(as.numeric(Plot_Area2),digits=2))

#-------------------#
# Estimate precision#
# pres/abs.         #
#-------------------#

RelAbund_precision <- data.frame(ranef(mdnc_lognormal_red)$Plot_Area2) %>% 
  rownames_to_column(.,var="Plot_Area2") %>%  
  mutate(Precision=1/Est.Error.Intercept,
         Plot_Area2=as.character(Plot_Area2)) %>%
  mutate(Plot_Area2=round(as.numeric(Plot_Area2),digits=2)) %>% 
  select(., Plot_Area2, Estimate.Intercept,Q2.5.Intercept,Q97.5.Intercept,Precision) %>% 
  left_join(.,plot_area_n,by="Plot_Area2")

#------------#
# Visualise  #
#------------#

PA_fun<-left_join(PA_accuracy,PA_precision,by="Plot_Area2")

# PA

est_PA<-ggplot(PA_fun, aes(x=Plot_Area2,y=Estimate.Intercept))+
  geom_hline(yintercept = 0,linetype=2)+
  geom_point(aes(size=log(N)))+
  geom_pointrange(aes(y=Estimate.Intercept, ymin=Q2.5.Intercept,ymax=Q97.5.Intercept))+
  labs(x="Plot area (m2)",y="Estimated intercept")+
  theme_bw()+theme(axis.title.y=element_text(colour="black",face="bold",size=8),
                   axis.title.x=element_text(colour="black",face="bold",size=8),
                   axis.text.x=element_text(colour="black",face="bold",size=8),
                   axis.text.y=element_text(colour=c("black"),face="bold",size=8),
                   plot.margin=margin(t=0.5,r=0.1,b=0.1,l=0.1,unit="cm"),
                   legend.position = "none",
                   panel.grid.major = element_blank(), panel.grid.minor = element_blank())

V_A_PA<-ggplot(PA_fun, aes(y=accurancy_xbar,x=Precision))+
  #geom_hline(yintercept = 0,linetype=2)+
  geom_point(aes(size=Plot_Area2))+
  geom_pointrange(aes(y=accurancy_xbar, ymin=accuracy_lwr,ymax=accuracy_upr))+
  labs(y="Accuracy of estimated probability of establishment",x="Precision (1/V)")+
  theme_bw()+theme(axis.title.y=element_text(colour="black",face="bold",size=8),
                   axis.title.x=element_text(colour="black",face="bold",size=8),
                   axis.text.x=element_text(colour="black",face="bold",size=8),
                   axis.text.y=element_text(colour=c("black"),face="bold",size=8),
                   plot.margin=margin(t=0.5,r=0.1,b=0.1,l=0.1,unit="cm"),
                   legend.position = "none",
                   panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# RelAbund

RelAbund_fun<-left_join(RelAbund_accuracy,RelAbund_precision,by="Plot_Area2")

est_RelAbund<-ggplot(RelAbund_fun, aes(x=Plot_Area2,y=Estimate.Intercept))+
  geom_hline(yintercept = 0,linetype=2)+
  geom_point(aes(size=log(N)))+
  geom_pointrange(aes(y=Estimate.Intercept, ymin=Q2.5.Intercept,ymax=Q97.5.Intercept))+
  labs(x="Plot area (m2)",y="Estimated intercept")+
  theme_bw()+theme(axis.title.y=element_text(colour="black",face="bold",size=8),
                   axis.title.x=element_text(colour="black",face="bold",size=8),
                   axis.text.x=element_text(colour="black",face="bold",size=8),
                   axis.text.y=element_text(colour=c("black"),face="bold",size=8),
                   plot.margin=margin(t=0.5,r=0.1,b=0.1,l=0.1,unit="cm"),
                   legend.position = "none",
                   panel.grid.major = element_blank(), panel.grid.minor = element_blank())

V_A_RelAbund<-ggplot(RelAbund_fun, aes(y=accurancy_xbar,x=Precision))+
  geom_hline(yintercept = 0,linetype=2)+
  geom_point(aes(size=Plot_Area2))+
  geom_pointrange(aes(y=accurancy_xbar, ymin=accuracy_lwr,ymax=accuracy_upr))+
  labs(y="Accuracy of estimated relative abundance",x="Precision (1/V)")+
  theme_bw()+theme(axis.title.y=element_text(colour="black",face="bold",size=8),
                   axis.title.x=element_text(colour="black",face="bold",size=8),
                   axis.text.x=element_text(colour="black",face="bold",size=8),
                   axis.text.y=element_text(colour=c("black"),face="bold",size=8),
                   plot.margin=margin(t=0.5,r=0.1,b=0.1,l=0.1,unit="cm"),
                   legend.position = "none",
                   panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# together
togg<-cowplot::plot_grid(est_PA,est_RelAbund,V_A_PA,V_A_RelAbund,
                         labels=c("a) Establishment (Species)","b) Relative abundance (Species)",
                                  "c) Establishment (Species)","d) Relative abundance (Species)"),
                         label_fontface = "bold", label_size = 7,
                         align="h", axis="tblr", ncol=2, nrow=2, greedy=TRUE)

ggsave(togg,filename = file.path("Figures", "PrecisionAccuracy_PA_Abund_Spp.png"), 
       width    = 20, 
       height   = 20, 
       units    = "cm",
       dpi=500)
