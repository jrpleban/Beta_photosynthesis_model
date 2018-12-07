## Script for setting up Bayesian Estimatation of Photosynthesis parameters from A/Ci & Light Responce curves w/ fluorescence parameter phiPS2
# lateest R version used with 3.5.1
# Updated: 10_19_2018. Jonathan R Pleban
### Script depends on  model for photosynthesis using combined Gas-exchange and Chlorophyll
##        fluorescence under limiting light conditions (light < 200 umol phonton m-2s-1)
## The model is based upon the  work of Yin et al (2009) PCE
## "Using combined measurements of gas exchange and
#      chlorophyll fluorescence to estimate parameters of a
#      biochemical C3 photosynthesis model: a critical appraisal
#      and a new integrated approach applied to leaves in a
#      wheat (Triticum aestivum) canopy"   and some other associate papers from Xinyou Yin
## Ccrit informed by data of Ci and Fv'Fm in external model ###
# Ccrit support based on Relialbe Estimation of biochemical paramters of C3 leaf...
###          Gu et al (2010) PCE Figs. 8,9,10
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
datalist1<-list(N=c(ACiN+LRN), N_ll=LRllN,
geno=c(ACi$Plant_grpF,LR$Plant_grpF),
geno_ll=LRll$Plant_grpF, Ngeno=Ngeno,
An=c(ACi$Photo,LR$Photo),
CiP=c(ACi$CP,LR$CP),
O=c(ACi$O,LR$O),
Q=c(ACi$PARi,LR$PARi),
T=c(ACi$Tleaf,LR$Tleaf),
A_ll=LRll$Photo,
phi2_ll=LRll$PhiPS2,
Inc_ll=LRll$PARi,
Kref=Kref, Tref=Tref, R=R, CCrit = 22 )
#############################
##### Parameter Set up ######
#############################
### parnames is used later to pull out parameters post model --- it is model specific as each has inherent complexity
parnames<-c("gm25", "Vcmax25", "gammaS25", "Rd25" ,
"phi2ll", "Jmax25", "s1" , "thetaJ", "Kc25", "Ko25","Egm", "EVcmax", "EgammaS", "ERd", "EKc", "EKo", "EJmax")
#### Set up for rjags #########
parameters = c(parnames,"mu.gm25", "mu.Vcmax25",  "mu.Rd25" ,
"mu.phi2ll", "mu.Jmax25",  "mu.s1","mu.thetaJ","tau.gm25", "tau.Vcmax25", "tau.Rd25" ,
"tau.phi2ll", "tau.Jmax25",  "tau.s1","tau.thetaJ","tau","tau_ll")
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
### This is the Yin model version
source("Mod_scripts/GECF_hier_model.R")
#################################
### running each curve
print("initialize models")
model1 <- jags.model(textConnection(GECF),
                     data = datalist1, n.chains=nChains , n.adapt=adaptSteps)
#################################
print("updating")
update(model1, burnInSteps) # Burnin for burnInSteps samples
##########################################
print("sampling chains")
##### mcmc_samples  #####
mcmc_samples1<- coda.samples(model1,
variable.names=parameters,
n.iter=nPerChain , thin=thinSteps )
### Convert mcmc_samples to a matrix
mcmcChain = as.matrix( mcmc_samples1)
# Convert precision (tau) to SD ###
sigma =1  / sqrt( mcmcChain[, "tau" ] )
hist(sigma)
mcmcChain = as.data.frame(cbind( mcmcChain, sigma ))
## calc gelman convergence diagnostic 
g1<-gelman.diag(mcmc_samples1,multivariate = FALSE)


Meds<-as.data.frame(apply(mcmcChain,2,median))
meds<-cbind(Meds,c(g1$psrf[,1],"NA"),c(g1$psrf[,2],"NA"))
colnames(meds)<-c("med","g","gmax")


print("writing samples")
setwd("~/Documents/Drought_Photosynthesis/Post_data")
write.table(mcmcChain,file=paste("mcmcChain_GECF_fixedCCrit_ACI&LR_no_late",Sys.Date(),sep = "_"), sep="\t", col.name=TRUE)
write.table(meds,file=paste("meds_GECF_fixedCCrit_ACI&LR_no_late",Sys.Date(),sep = "_"), sep="\t", col.name=TRUE)


