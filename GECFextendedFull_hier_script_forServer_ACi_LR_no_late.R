## Bayesian Estimate of PS paramters from Brassica ACi & LR curves ###
 #w/ fluorescence paramters Fv'Fm', and phiPS2
## Ccrit informed by data of Ci and Fv'Fm ###
# R 3.2.1
# Updated: 11_20_2018. Jonathan R Plean.
########################################################################
### Script depends on  model for photosyntheis using combined
##    Gas-exchange and chlorophyll fluorescence to inform description of electron transport rate
##    making use of PhiPS2 vs light (Q) relationship across Q from 0-2000 umol photon m-2 s-1
###   PhiPS2 vs Q modeled as an expodential decay fucntion where betaPSII is the decay rate as light approaches infinity
###
##  rjags model conded as GECfextendedFull_hier_model_Ccrit
#   Model is largely based on Yin 2009 PCE
# Ccrit support based on Relialbe Estimation of biochemical paramters of C3 leaf...
###          Gu et al (2010) PCE Figs. 8,9,10
### use of Phi2 vs PAR relationship developed by JRPleban for this work
########################################################################
########################################################################
library("rjags") ### package for Gibbs Sampling (JAGS)
#### working directoy for brassica drought project July-Aug 2017 data growth Chamber
setwd("~/Documents/Drought_Photosynthesis")
######################
#####DATA SET UP ######
F<-read.delim("data/Full_data_compiled_Pleban_Guadagno_Summer2017.txt")
###
## seperate LR and ACi curves
LR1<-F[F$curve=="LR",]
## remove late drought as no Associated ACi data
LR<-LR1[!LR1$drought_L=="late",]
ACi<-F[F$curve=="Aci",]
### low light data only
LRll<-LR[LR$PARi<=205,]
### constants
Kref=298.15; R=0.008314; Tref=25
## hier input needs
ACiN<-length(ACi[,1])
LRN<-length(LR[,1])
Ngeno=12   # this is the number of genotypes (4) times the number of treatments (3)
#############################
### pass data into models as a datalist with hierarchy described by Plant_grpF
datalist1<-list(N=ACiN+LRN,N_L=LRN,Ngeno=Ngeno,
geno=c(ACi$Plant_grpF, LR$Plant_grpF),
geno_ll=LR$Plant_grpF,
An=c(ACi$Photo,LR$Photo),
CiP=c(ACi$CP,LR$CP),
O=c(ACi$O,LR$O),
Q=c(ACi$PARi,LR$PARi),
T=c(ACi$Tleaf,LR$Tleaf),
A_L=LR$Photo,
phi2=LR$PhiPS2,
Inc_L=LR$PARi,
Kref=Kref, Tref=Tref, R=R, CCrit=22 )
#############################
##### Parameter Set up ######
#############################
### parnames is used later to pull out parameters post model --- it is model specific as each has inherent complexity
parnames<-c("gm25", "Vcmax25", "gammaS25", "Rd25" ,
            "Kc25", "Ko25",
            "Egm", "EVcmax", "EgammaS", "ERd", "EKc", "EKo", "EJmax",
            "s1" ,"alpha", "beta", "kappa")
#### Set up for rjags #########
parameters = c(parnames,"mu.gm25", "mu.Vcmax25",  "mu.Rd25" ,
  "mu.s1","mu.alpha", 
"mu.beta","mu.kappa","tau.gm25",
"tau.Vcmax25",  "tau.Rd25" ,
  "tau.CCrit", "tau.s1",
"tau.kappa","tau.alpha","tau.beta","tau","tau_ll")
#############################
#####  MCMC  Set up    ######
#############################
### chains and burnin
adaptSteps = 10000             # Number of steps to "tune" the samplers.
burnInSteps = 50000            # Number of steps to "burn-in" the samplers.
nChains = 4                   # Number of chains to run.
DICsteps= 20000                # Number of steps of sample DIC
numSavedSteps= 2500        # Total number of steps in chains to save.
thinSteps=20                   # Number of steps to "thin" (1=keep every step).
nPerChain = ceiling( ( numSavedSteps * thinSteps ) / nChains ) # Steps per chain.
##################################
##################################
##### Impliment model in JAGS ####
source("Mod_scripts/GECFextendedfull_hier_model.R")
################################
### running each curve
print("initialize models")
model1 <- jags.model(textConnection(GECFext), inits = inits,
                     data = datalist1, n.chains=nChains , n.adapt=adaptSteps)
#################################
print("updating")
update(model1, burnInSteps) # Burnin for burnInSteps samples
#################################
print("sampling chains")
##### mcmc_samples  #####
mcmc_samples1<- coda.samples(model1,
variable.names=parameters,
n.iter=nPerChain , thin=thinSteps )
####### Plot results #####
#plot(mcmc_samples1)
mcmcChain = as.matrix( mcmc_samples1)
chainLength = NROW(mcmcChain)
# Convert precision (tau) to SD###
sigma =1  / sqrt( mcmcChain[, "tau" ] )
#hist(sigma)
mcmcChain = as.data.frame(cbind( mcmcChain, sigma ))
g1<-gelman.diag(mcmc_samples1,multivariate = FALSE)
##########################################
#######  median estimates ############
Meds<-as.data.frame(apply(mcmcChain,2,median))
meds<-cbind(Meds,c(g1$psrf[,1],"NA"),c(g1$psrf[,2],"NA"))
colnames(meds)<-c("med","g","gmax")
##########################################
####### Save posteriors ############
print("writing samples")
setwd("~/Documents/Drought_Photosynthesis/Post_data")
write.table(mcmcChain,file=paste("mcmcChain_Fullyextended_fixedCCrit_ACI&LR_no_late",Sys.Date(),sep = "_"), sep="\t", col.name=TRUE)
write.table(meds,file=paste("meds_Fullyextended_fixedCCrit_ACI&LR_no_late",Sys.Date(),sep = "_"), sep="\t", col.name=TRUE)


