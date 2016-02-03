SIH_function<-function(dispersal=0.001,model_ = NA,species=9,patches=30,eff_vary=F,eff_sd=NA, eff_mean = NA){ 
	#can define things here and then don't need to define them below 
  require(cluster)
  require(picante)
  
  #Constants####
  N<- matrix(10,ncol=species,nrow=patches) # Community x Species abundance matrix
  R<-rep(10*(species/10),patches) #Initial resources
  N0<-N
  R0<-R
  
  rInput<-150 #resource input
  rLoss<-10 #resource loss 
  mort<-0.2 #mortality
  Ext<- 0.1 #extinction Threshold
  eff<-0.2 #conversion efficiency
  
  ePeriod<-40000 #period of env sinusoidal fluctuations
  eAMP<-1 #amplitude of envrionment sinusoidal fluctuations
  
  Tmax<-140000 #number of time steps in Simulation
  DT<- 0.08 # % size of discrete "time steps" - this is the Euler value
  
  #simulating phylogeny
  uncorr <- matrix(data = c(1,0,0,0,1,0,0,0,1), nrow = 3, ncol = 3)
  phylo <- rcoal(species)
  
  #phylo<-as.phylo(hclust(daisy(cbind(eOptimum,eff_values)),method ="ward"))
  #plot(phylo,show.tip.label = F)
  #tiplabels(pch=22,bg=heat.colors(species)[1:species]) #heat.colors = white -> red

if(model_ == "BM"){
#getting trait values
traits<-sim.char(phylo,uncorr,nsim = 1, model = "BM", root = 1) #traits is an array, oddly enough
traits.stand<-decostand(traits[,1,1],"range") #normalizes the H trait b/w 0 and 1
traits.stand2<-decostand(traits[,2,1],"range") #normalizes the second H trait between 0 and 1
traits.stand3<-decostand(traits[,3,1],"standardize")
traits.stand3_fixed <- ((eff_sd)*traits.stand3)+(eff_mean) #makes e a normal r.v. with mean = 0.2, s.d. = 0.005
traitsfinal<-cbind(traits.stand, traits.stand2, traits.stand3_fixed)
}  

if(model_ == "random"){
#getting trait values
traits<-sim.char(phylo,uncorr,nsim = 1, model = "BM", root = 1) #traits is an array, oddly enough
traits.stand<-decostand(traits[,1,1],"range") #normalizes the H trait b/w 0 and 1
traits.stand2<-decostand(traits[,2,1],"range") #normalizes the second H trait between 0 and 1
traits.stand3<-decostand(traits[,3,1],"standardize")
traits.stand3_fixed <- ((eff_sd)*traits.stand2)+(eff_mean) #makes e a normal r.v. with mean = 0.2, s.d. = 0.005
traitsfinal<-cbind(traits.stand, traits.stand2, traits.stand3_fixed)
traitsfinal<-traitsfinal[sample(nrow(traitsfinal)),]

} 
###NEED TO FIX BELOW HERE STILL###
#trait values -> model parameters
#species environmental optima (H)
if(eff_vary==T){eOptimum1<-traitsfinal[,1]
eOptimum2<-traitsfinal[,2]
} else{
	eOptimum1<-1-seq(0,eAMP, by=eAMP/(species-1))
	eOptimum2<-1-seq(0,eAMP, by=eAMP/(species-1)) #this might not work
	}

#species consumption efficiency - modified!!!	
#species consumption efficiency (below)
eff<-rep(eff,species) #still a vector, but all just one value now 
  

  #dispersal conditions####
  dispersal_matrix<-matrix(1/(patches-1),patches,patches)
  diag(dispersal_matrix)<-0
  
  calc.immigration <- function(N,a,dispersal_matrix) dispersal_matrix%*%N*rep(a,each=patches)
  
  Prod<-matrix(NA,species*patches,40000)
  Abund<-Prod
  
   Meta_dyn<-data.frame(Species_sorting=rep(NA,40000),Mass_effects=NA,Base_growth=NA)
  Species_data<-array(NA,dim=c(40000,species,2),dimnames = list(1:40000,1:species,c("Abundance","Occupancy")))
  
  for(TS in 1:Tmax){
    Immigrants<-calc.immigration(N,dispersal,dispersal_matrix)
    envt.v<-0.5*eAMP*(sin((2*pi/ePeriod)*TS+1+(1:patches)*2*pi/patches)+1)
    consume <- 0.05*((1.5-abs(sapply(eOptimum1,'-',envt.v))) + (1.5-abs(sapply(eOptimum2,'-',envt.v))))
    Nt <- N*(1+DT*(rep(eff,each=patches)*R*consume - dispersal - mort)) + DT*Immigrants #abundance step
    Immigrants0<-calc.immigration(N0,0,dispersal_matrix)
    Nt0 <- N0*(1+DT*(rep(eff,each=patches)*R0*consume - 0 - mort)) + DT*Immigrants0

    Rt <- DT*rInput+R*(1-DT*(rLoss + rowSums(consume*N))) #resource step 
    Rt0 <- DT*rInput+R0*(1-DT*(rLoss + rowSums(consume*N0))) #resource step   
    N <- Nt * (Nt>Ext) # set to 0 if below extinction threshold
    R <- Rt
    
    if(TS>=100000){
      Prod[,(TS-100000)] <- c(t(eff*consume*R*N))
      Abund[,(TS-100000)] <- c(t(N))
      
      fitness<-((N*(1+DT*(rep(eff,each=patches)*R*consume - dispersal - mort)))-N)*(Nt>Ext)
      fitness_w_disp<-((N*(1+DT*(rep(eff,each=patches)*R*consume - dispersal - mort)) + DT*Immigrants)-N)*(Nt>Ext)
      fitness0<-(N0*(1+DT*(rep(eff,each=patches)*R0*consume - mort))-N0)*(Nt0>Ext)
      home_prod<-mean(rowSums(fitness_w_disp*(fitness>0)))
      disp_prod_ME<-mean(rowSums(fitness_w_disp*(fitness<0 & fitness_w_disp>=0)))
      
      base_prod<-mean(rowSums(fitness0*(fitness0>0)))
      total_prod<-home_prod+disp_prod_ME
      
      home_prod_prop<-home_prod/total_prod
      SS_prod<-home_prod-base_prod
      SS_prod[SS_prod<0]<-0
      if(mean(rowSums(N>0))<=1){SS_prod<-0}
      SS<-(SS_prod/home_prod)*home_prod_prop
      SS[is.nan(SS)]<-0
      if(total_prod==0){SS<-NA}
      Meta_dyn$Species_sorting[(TS-100000)]<-SS
      
      ME<-(disp_prod_ME)/total_prod
      ME[is.nan(ME)]<-0
      if(total_prod==0){ME<-NA}
      Meta_dyn$Mass_effects[(TS-100000)]<-ME
      
      BP<-home_prod_prop*(1-(SS_prod/home_prod))
      BP[is.nan(BP)]<-0
      if(total_prod==0){BP<-NA}
      Meta_dyn$Base_growth[(TS-100000)]<-BP
      
      Species_data[(TS-100000),,1]<-colSums(N)
      Species_data[(TS-100000),,2]<-colSums(N>0)
    }
    N <- Nt * (Nt>Ext) # set to 0 if below extinction threshold
    R <- Rt
    
    N0 <- Nt0 * (Nt0>Ext) # set to 0 if below extinction threshold
    R0 <- Rt0
  } 
  
  Prod<-array(t(Prod),dim=c(40000,species,patches))
  Prod<-Prod[seq(1,40000,100),,]
  
  Abund<-array(t(Abund),dim=c(40000,species,patches))
  Abund<-Abund[seq(1,40000,100),,] #this only takes a subset, only every 100 time steps 
  return(list(Prod=Prod,Abund=Abund,phylo=phylo,Meta_dyn=Meta_dyn,Spec_data=apply(Species_data,3,colMeans))) #will only return these things, nothing else
}


