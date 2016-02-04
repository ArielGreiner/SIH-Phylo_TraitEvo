require(dplyr) 
Data_storage_total<-summarise(group_by(Data_storage, Dispersal, TEvoModel, Scale), Mean_SR=mean(SR,na.rm=T), SD_SR=sd(SR,na.rm=T), Mean_Biomass=mean(Biomass,na.rm=T), SD_Biomass=sd(Biomass,na.rm=T), Mean_PD=mean(PD,na.rm=T), SD_PD=sd(PD,na.rm=T), Mean_MPD_abund=mean(MPD_abund,na.rm=T), SD_MPD_abund=sd(MPD_abund,na.rm=T), Mean_MPD_pa=mean(MPD_pa,na.rm=T), SD_MPD_pa=sd(MPD_pa,na.rm=T), Mean_MNTD_abund=mean(MNTD_abund, na.rm=T), SD_MNTD_abund=sd(MNTD_abund,na.rm=T), Mean_MNTD_pa=mean(MNTD_pa,na.rm=T), SD_MNTD_pa=sd(MNTD_pa,na.rm=T), Mean_Kstat=mean(Kstat, na.rm=T), SD_Kstat=sd(Kstat, na.rm=T), Mean_shannonhill=mean(shannonhill, na.rm=T), SD_shannonhill=sd(shannonhill, na.rm=T), Mean_shannonhillbeta=mean(shannonhillbeta, na.rm=T), SD_shannonhillbeta=sd(shannonhillbeta, na.rm=T), Mean_speciessorting=mean(speciessorting, na.rm=T), SD_speciessorting=sd(speciessorting, na.rm=T), Mean_masseffects=mean(masseffects, na.rm=T), SD_masseffects=sd(masseffects, na.rm=T), Mean_basegrowth = mean(basegrowth, na.rm=T), SD_basegrowth = sd(basegrowth, na.rm=T), Mean_sesMPD_abund_z = mean(sesMPD_abund_z,na.rm=T), SD_sesMPD_abund_z = sd(sesMPD_abund_z,na.rm=T),
Mean_sesMPD_abund_p = mean(sesMPD_abund_p,na.rm=T), SD_sesMPD_abund_p = sd(sesMPD_abund_p,na.rm=T), Mean_sesMNTD_abund_z = mean(sesMNTD_abund_z,na.rm=T), SD_sesMNTD_abund_z = sd(sesMNTD_abund_z,na.rm=T), Mean_sesMNTD_abund_p = mean(sesMNTD_abund_p,na.rm=T), SD_sesMNTD_abund_p = sd(sesMNTD_abund_p,na.rm=T), sum_phylogeven_mpd = sum(phylogeven_mpd), sum_phylogeven_mntd = sum(phylogeven_mntd), sum_phylogcluster_mpd = sum(phylogcluster_mpd), sum_phylogcluster_mntd = sum(phylogcluster_mntd), sum_pe_pc_mpd = sum(phylogeven_mpd + phylogcluster_mpd), sum_pe_pc_mntd = sum(phylogeven_mntd + phylogcluster_mntd))

require(ggplot2)
#MPD abund vs dispersal rates and then look at 3 different levels of relationship between traits and the phylogeny
ggplot(Data_storage_total,aes(x=Dispersal,y=Mean_SR,color=interaction(Scale, TEvoModel),group=interaction(Scale, TEvoModel),fill=interaction(Scale, TEvoModel),alpha=0.1))+
  geom_line(size=2)+ #plots data as lines
  geom_ribbon(aes(ymin=Mean_SR-SD_SR,ymax=Mean_SR+SD_SR),width=0.1)+
  scale_x_log10(breaks=DispV)+ #sets x axis to log10 scale
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines

#MPDabund vs species-sorting, MPDabund vs mass effects, MPDabund vs base growth, 
ggplot(Data_storage_total,aes(x=Mean_speciessorting,y=Mean_MPD_abund,group=TEvoModel,color=TEvoModel, fill=TEvoModel,alpha=0.1))+
  geom_line(size=1.5)+
  geom_ribbon(aes(ymin=Mean_MPD_abund-SD_MPD_abund,ymax=Mean_MPD_abund+SD_MPD_abund),width=0.1)+
  theme_bw(base_size = 15)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_x_log10(breaks=DispV)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


ggplot(Data_storage_total,aes(x=Dispersal,y=Mean_MPD_abund,group=interaction(TEvoModel, Scale),color=interaction(TEvoModel,Scale),fill=interaction(TEvoModel, Scale),alpha=0.1))+
  geom_line(size=1.5)+
geom_ribbon(aes(ymin=Mean_MPD_abund-SD_MPD_abund,ymax=Mean_MPD_abund+SD_MPD_abund),width=0.1)+
theme_bw(base_size = 15)+
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
scale_x_log10(breaks=DispV)+
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(Data_storage_total,aes(x=Dispersal,y=Mean_shannonhillbeta,group=interaction(TEvoModel, Scale),color=interaction(TEvoModel,Scale),fill=interaction(TEvoModel, Scale),alpha=0.1))+
  geom_line(size=1.5)+
geom_ribbon(aes(ymin=Mean_shannonhillbeta-SD_shannonhillbeta,ymax=Mean_shannonhillbeta + SD_shannonhillbeta),width=0.1)+
theme_bw(base_size = 15)+
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
scale_x_log10(breaks=DispV)+
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(Data_storage_total,aes(x=Dispersal,y=Mean_shannonhill,group= Scale,color=Scale,fill=Scale,alpha=0.1))+
  geom_line(size=1.5)+
geom_ribbon(aes(ymin=Mean_shannonhill-SD_shannonhill,ymax=Mean_shannonhill + SD_shannonhill),width=0.1)+
theme_bw(base_size = 15)+
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
scale_x_log10(breaks=DispV)+
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


ggplot(Data_storage_total,aes(x=Mean_SR,y=Mean_Biomass,color=factor(Dispersal),group=interaction(Scale, TEvoModel)))+
geom_point()+ #plots data as points
geom_path()+
geom_errorbar(aes(ymin=Mean_Biomass-SD_Biomass,ymax=Mean_Biomass+SD_Biomass),width=0.1)+
facet_grid(Scale~TEvoModel,scale="free")+
theme_bw(base_size = 18)+ #gets rid of grey background
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(Data_storage_total,aes(x=Mean_speciessorting,y=Mean_MPD_abund,color=factor(Dispersal),group=interaction(Scale, TEvoModel)))+

geom_point()+ #plots data as points
geom_path()+
geom_errorbar(aes(ymin=Mean_MPD_abund-SD_MPD_abund,ymax=Mean_MPD_abund+SD_MPD_abund),width=0.1)+

facet_grid(TEvoModel~.,scale="free")+

  theme_bw(base_size = 15)+

  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+

  scale_x_log10(breaks=DispV)+

  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  #Plotting p-value things relating to ses measures
ggplot(Data_storage_total,aes(x=Dispersal,y=Mean_sesMPD_abund_z,color=factor(sum_pe_pc_mpd)/100,group=interaction(Scale, TEvoModel)))+
  geom_line(size=1.5)+
geom_errorbar(aes(ymin=Mean_sesMPD_abund_z-SD_sesMPD_abund_z,ymax=Mean_sesMPD_abund_z + SD_sesMPD_abund_z),width=0.1)+
theme_bw(base_size = 15)+
facet_grid(Scale~TEvoModel,scale="free")+
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
scale_x_log10(breaks=DispV)+
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
 

#ignore under here, it doesn't work as it should...
Meta_dynamics<-data.frame(Mean_Proportion_of_Production=c(Data_storage_total$Mean_basegrowth,Data_storage_total$Mean_speciessorting,Data_storage_total$Mean_masseffects),SD_Proportion_of_Production=c(Data_storage_total$SD_basegrowth,Data_storage_total$SD_speciessorting,Data_storage_total$SD_masseffects),Dynamic=rep(c("Base growth","Species sorting","Mass effects"), each=length(DispV)*length(TraitEvo)), TraitEvoModel=rep(c("BM","random","conserved"), each=length(DispV)), Mean_MPD_abund = Data_storage_total$Mean_MPD_abund, SD_MPD_abund = Data_storage_total$SD_MPD_abund, Dispersal=rep(DispV,each=nreplicates*length(TraitEvo)),TEvoModel=factor(TraitEvo), ReplicateNum=rep(reps,each=length(TraitEvo)),
Scale=rep(c("Local","Regional"),each=length(DispV)*nreplicates*length(TraitEvo)))