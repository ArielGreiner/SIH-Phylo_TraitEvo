require(geiger)
nspecies<-7 #the number of species
npatches<-10 #the number of patches
nreplicates<-100 #number of replicates
nfunctions<-1

DispV<-c(0.0001,0.0005,0.001,0.005,0.01,0.05,0.1,0.5,1) #the dispersal rates
TraitEvo<-c("BM", "random")  
reps <- seq(from = 1, to = nreplicates, by = 1)

Data_storage<-data.frame(SR=NA,Biomass=NA,Biomass_CV=NA,PD=NA,MPD_abund=NA,MPD_pa=NA,MNTD_pa = NA, MNTD_abund=NA, Kstat = NA, shannonhill = NA, shannonhillbeta = NA, speciessorting = NA, masseffects = NA, basegrowth = NA, Dispersal=rep(DispV,each=nreplicates*length(TraitEvo)),TEvoModel=factor(TraitEvo), ReplicateNum=rep(reps,each=length(TraitEvo)),
Scale=rep(c("Local","Regional"),each=length(DispV)*nreplicates*length(TraitEvo))) #building the data frame

MTraits<-t(matrix(1,nspecies,nfunctions))
for(l in 1:length(TraitEvo)){
	if(l == 2){
		print(l)
	}
for(j in 1:nreplicates){
		if(j == 25){
		print(j)
	}
	if(j == 50){
		print(j)
	}
	
		if(j == 25){
		print(j)
	}
#runs the SIH model at all dispersal rates in DispV and saves the abundances and productivity in a list
set.seed(j)
eff_sd <- 0.005
eff_mean <- 0.2
SIH_data<-sapply(DispV,model_ = TraitEvo[l], SIH_function,species=nspecies,patches=npatches,eff_vary=T,eff_sd=eff_sd, eff_mean = eff_mean)

shannon_alpha <- matrix(data = rep(0,npatches*length(DispV)), nrow = npatches, ncol = length(DispV))
shannon_gamma <- rep(0, length = nspecies)
for(i in 1:length(DispV)){
 #BIOMASS
  LFunc_rate1<-apply(SIH_data[["Abund",i]],1,function(M){MTraits%*%M}) #calculates the amount of each function in each patch based on the abundances from the SIH model
  LFunc_rate<-array(t(LFunc_rate1),dim=c(nrow(SIH_data[["Abund",1]]),nfunctions,npatches)) #arranges the functions into an array that is easier to use
  Data_storage$Biomass[Data_storage$Dispersal==DispV[i] & Data_storage$TEvoModel==TraitEvo[l] & Data_storage$ReplicateNum==j &  Data_storage$Scale == "Local"]<-mean(apply(LFunc_rate,3,colMeans))
  Data_storage$Biomass[Data_storage$Dispersal==DispV[i] & Data_storage$TEvoModel==TraitEvo[l] & Data_storage$ReplicateNum==j & Data_storage$Scale == "Regional"]<-colMeans(apply(LFunc_rate,2,rowSums)) #calculates and saves the mean regional rate for each function
  
	#calculate species richness at the local and regional scale	
  Data_storage$SR[Data_storage$Dispersal==DispV[i] & Data_storage$TEvoModel==TraitEvo[l] & Data_storage$ReplicateNum==j &  Data_storage$Scale == "Local"] <-mean(rowMeans(apply(SIH_data[["Abund",i]]>0,3,rowSums))) 
  Data_storage$SR[Data_storage$Dispersal==DispV[i] & Data_storage$TEvoModel==TraitEvo[l] & Data_storage$ReplicateNum==j & Data_storage$Scale == "Regional"]<-mean(rowSums(apply(SIH_data[["Abund",i]],2,rowMeans)>0)) #regional SR of all species in each time step
  
  #At the local scale...
    interspec_mat <- cophenetic(SIH_data[["phylo",i]])
     colnames(interspec_mat) <- 1:nspecies
     rownames(interspec_mat) <- 1:nspecies	
  com_data<-t(SIH_data[["Abund",i]][400,,])
  colnames(com_data)<-1:nspecies
  SIH_data[["phylo",i]]$tip.label <- 1:nspecies
  x = interspec_mat
  phy = SIH_data[["phylo", i]]
  Data_storage$Kstat[Data_storage$Dispersal==DispV[i] & Data_storage$TEvoModel==TraitEvo[l] & Data_storage$ReplicateNum==j & Data_storage$Scale == "Regional"]<- mean(Kcalc(x=interspec_mat,phy=SIH_data[["phylo",i]])) #calculating K statistic, not sure if should take the mean of all of the values 
  
  Data_storage$PD[Data_storage$Dispersal==DispV[i] & Data_storage$TEvoModel==TraitEvo[l] & Data_storage$ReplicateNum==j & Data_storage$Scale == "Local"]<-mean(pd(com_data,SIH_data[["phylo",i]])$PD) 
  Data_storage$MPD_abund[Data_storage$Dispersal==DispV[i] & Data_storage$TEvoModel==TraitEvo[l] & Data_storage$ReplicateNum==j & Data_storage$Scale == "Local"]<-mean(mpd(com_data,cophenetic(SIH_data[["phylo",i]]),abundance.weighted = T)) 
  Data_storage$MPD_pa[Data_storage$Dispersal==DispV[i] & Data_storage$TEvoModel==TraitEvo[l] & Data_storage$ReplicateNum==j & Data_storage$Scale == "Local"]<-mean(mpd(com_data,cophenetic(SIH_data[["phylo",i]]),abundance.weighted = F))
  Data_storage$MNTD_abund[Data_storage$Dispersal==DispV[i] & Data_storage$TEvoModel==TraitEvo[l] & Data_storage$ReplicateNum==j & Data_storage$Scale == "Local"]<-mean(mntd(com_data,cophenetic(SIH_data[["phylo",i]]),abundance.weighted = T))
   Data_storage$MNTD_pa[Data_storage$Dispersal==DispV[i] & Data_storage$TEvoModel==TraitEvo[l] & Data_storage$ReplicateNum==j & Data_storage$Scale == "Local"]<-mean(mntd(com_data,cophenetic(SIH_data[["phylo",i]]),abundance.weighted = F))
   
#alpha, beta, gamma
for(k in 1:npatches){
	for(m in 1:nspecies){
		relabund <- com_data[k,m]/sum(com_data[k,])
		shannon_alpha[k,i][is.na(shannon_alpha[k,i])] <- 0 #checks if NA and replaces with 0
		shannon_alpha[k,i] <- relabund*log(relabund) + shannon_alpha[k,i]
		gammaabund <- sum(com_data[,m])/sum(com_data)
		shannon_gamma[i][is.na(shannon_gamma[i])] <- 0
		shannon_gamma[i] <- gammaabund*log(gammaabund) + shannon_gamma[i]
	}	
}
avgshannon_a <- prod(exp(shannon_alpha[,i]))^(1/length(shannon_alpha[,i]))
Data_storage$shannonhill[Data_storage$Dispersal==DispV[i] & Data_storage$TEvoModel==TraitEvo[l] & Data_storage$ReplicateNum==j & Data_storage$Scale == "Local"] <- avgshannon_a	
Data_storage$shannonhill[Data_storage$Dispersal==DispV[i] & Data_storage$TEvoModel==TraitEvo[l] & Data_storage$ReplicateNum==j & Data_storage$Scale == "Regional"] <- exp(shannon_gamma[i])
Data_storage$shannonhillbeta[Data_storage$Dispersal==DispV[i] & Data_storage$TEvoModel==TraitEvo[l] & Data_storage$ReplicateNum==j & Data_storage$Scale == "Regional"]<- exp(shannon_gamma[i])/avgshannon_a
   
    #At the regional scale
  com_data<-matrix(colSums(t(SIH_data[["Abund",i]][400,,])),1,nspecies)
  colnames(com_data)<-1:nspecies
  Data_storage$PD[Data_storage$Dispersal==DispV[i] & Data_storage$TEvoModel==TraitEvo[l] & Data_storage$ReplicateNum==j & Data_storage$Scale == "Regional"]<-pd(com_data,SIH_data[["phylo",i]])$PD
  Data_storage$MPD_abund[Data_storage$Dispersal==DispV[i] & Data_storage$TEvoModel==TraitEvo[l] & Data_storage$ReplicateNum==j & Data_storage$Scale == "Regional"]<-mpd(com_data,cophenetic(SIH_data[["phylo",i]]),abundance.weighted = T)
  Data_storage$MPD_pa[Data_storage$Dispersal==DispV[i] & Data_storage$TEvoModel==TraitEvo[l] & Data_storage$ReplicateNum==j & Data_storage$Scale == "Regional"]<-mpd(com_data,cophenetic(SIH_data[["phylo",i]]),abundance.weighted = F)
  Data_storage$MNTD_abund[Data_storage$Dispersal==DispV[i] & Data_storage$TEvoModel==TraitEvo[l] & Data_storage$ReplicateNum==j & Data_storage$Scale == "Regional"]<-mntd(com_data,cophenetic(SIH_data[["phylo",i]]),abundance.weighted = T)  
   Data_storage$MNTD_pa[Data_storage$Dispersal==DispV[i] & Data_storage$TEvoModel==TraitEvo[l] & Data_storage$ReplicateNum==j & Data_storage$Scale == "Regional"]<-mntd(com_data,cophenetic(SIH_data[["phylo",i]]),abundance.weighted = F)
   
   #ss, me, bg
   Data_storage$speciessorting[Data_storage$Dispersal==DispV[i] & Data_storage$TEvoModel==TraitEvo[l] & Data_storage$ReplicateNum==j & Data_storage$Scale == "Regional"]<-SIH_data[["Meta_dyn",i]][400,,]$Species_sorting
    Data_storage$masseffects[Data_storage$Dispersal==DispV[i] & Data_storage$TEvoModel==TraitEvo[l] & Data_storage$ReplicateNum==j & Data_storage$Scale == "Regional"]<-SIH_data[["Meta_dyn",i]][400,,]$Mass_effects
      Data_storage$basegrowth[Data_storage$Dispersal==DispV[i] & Data_storage$TEvoModel==TraitEvo[l] & Data_storage$ReplicateNum==j & Data_storage$Scale == "Regional"]<-SIH_data[["Meta_dyn",i]][400,,]$Base_growth 
    
}
}
}

#going to have to fix the total data storage thing big time
