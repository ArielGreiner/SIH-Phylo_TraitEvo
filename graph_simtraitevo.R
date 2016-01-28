require(dplyr) 
Data_storage_total<-summarise(group_by(Data_storage, Dispersal, TEvoModel, Scale), Mean_SR=mean(SR,na.rm=T), SD_SR=sd(SR,na.rm=T), Mean_Biomass=mean(Biomass,na.rm=T), SD_Biomass=sd(Biomass,na.rm=T), Mean_PD=mean(PD,na.rm=T), SD_PD=sd(PD,na.rm=T), Mean_MPD_abund=mean(MPD_abund,na.rm=T), SD_MPD_abund=sd(MPD_abund,na.rm=T), Mean_MPD_pa=mean(MPD_pa,na.rm=T), SD_MPD_pa=sd(MPD_pa,na.rm=T), Mean_MNTD_abund=mean(MNTD_abund, na.rm=T), SD_MNTD_abund=sd(MNTD_abund,na.rm=T), Mean_MNTD_pa=mean(MNTD_pa,na.rm=T), SD_MNTD_pa=sd(MNTD_pa,na.rm=T), Mean_Kstat=mean(Kstat, na.rm=T), SD_Kstat=sd(Kstat, na.rm=T), Mean_shannonhill=mean(shannonhill, na.rm=T), SD_shannonhill=sd(shannonhill, na.rm=T), Mean_shannonhillbeta=mean(shannonhillbeta, na.rm=T), SD_shannonhillbeta=sd(shannonhillbeta, na.rm=T), Mean_speciessorting=mean(speciessorting, na.rm=T), SD_speciessorting=sd(speciessorting, na.rm=T), Mean_masseffects=mean(masseffects, na.rm=T), SD_masseffects=sd(masseffects, na.rm=T), Mean_basegrowth = mean(basegrowth, na.rm=T), SD_basegrowth = sd(basegrowth, na.rm=T))

#MPD abund vs dispersal rates and then look at 3 different levels of relationship between traits and the phylogeny
ggplot(Data_storage_total,aes(x=Dispersal,y=Mean_MPD_abund,group=TEvoModel,color=TEvoModel, fill=TEvoModel,alpha=0.1))+
  geom_line(size=1.5)+
  geom_ribbon(aes(ymin=Mean_MPD_abund-SD_MPD_abund,ymax=Mean_MPD_abund+SD_MPD_abund),width=0.1)+
  theme_bw(base_size = 15)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_x_log10(breaks=DispV)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#MPDabund vs species-sorting, MPDabund vs mass effects, MPDabund vs base growth, 
ggplot(Data_storage_total,aes(x=Mean_speciessorting,y=Mean_MPD_abund,group=TEvoModel,color=TEvoModel, fill=TEvoModel,alpha=0.1))+
  geom_line(size=1.5)+
  geom_ribbon(aes(ymin=Mean_MPD_abund-SD_MPD_abund,ymax=Mean_MPD_abund+SD_MPD_abund),width=0.1)+
  theme_bw(base_size = 15)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_x_log10(breaks=DispV)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
 dev.off() 

#ignore under here
Meta_dynamics<-data.frame(Mean_Proportion_of_Production=c(Data_storage_total$Mean_basegrowth,Data_storage_total$Mean_speciessorting,Data_storage_total$Mean_masseffects),SD_Proportion_of_Production=c(Data_storage_total$SD_basegrowth,Data_storage_total$SD_speciessorting,Data_storage_total$SD_masseffects),Dynamic=rep(c("Base growth","Species sorting","Mass effects"), each=length(DispV)*length(TraitEvo)), TraitEvoModel=rep(c("BM","random","conserved"), each=length(DispV)), Mean_MPD_abund = Data_storage_total$Mean_MPD_abund, SD_MPD_abund = Data_storage_total$SD_MPD_abund, Dispersal=rep(DispV,each=nreplicates*length(TraitEvo)),TEvoModel=factor(TraitEvo), ReplicateNum=rep(reps,each=length(TraitEvo)),
Scale=rep(c("Local","Regional"),each=length(DispV)*nreplicates*length(TraitEvo)))