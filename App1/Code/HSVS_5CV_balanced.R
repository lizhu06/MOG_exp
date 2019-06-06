########################
# SOG predict ER
# 5 cross-validation
# 12/09/2016
#######################
rm(list=ls())
library(glmnet)
library(grpreg)
library(SGL)
library(AUC)
library(RcppEigen)
library(Rcpp)
library(snowfall)

##### duplicate genes belonging to multiple groups
setwd("/mnt/glusterfs/liz86/MOG_Regression/Exp/app1_8pathways/")
load("Data/MOG_ER_CV_Pre_balanced.RData")  # only use fold-id
load("Data/SOG_ER_CV_Pre_balanced.RData")
#nsim <- 2000
#start_nsim <- 1001

## function for cross-validation

HSVS_CV <- function(X_s, Y, g_index, 
	  types, nsim, start_nsim, seed, fold_id, fold) {
	source("/mnt/glusterfs/liz86/MOG_Regression/R/HSVS/HSVS_modify_v3_binary.R")
	X_s_train <- X_s[which(fold_id != fold), ]
	Y_train <- Y[which(fold_id != fold)]
	res <- HSVS_binary(burnin=(start_nsim-1),N=(nsim-start_nsim+1),Y_train,X_s_train, gindex=g_index)  # N is the number of ite after burn-in to be saved
	save(res, file=paste("/mnt/glusterfs/liz86/MOG_Regression/Exp/app1_8pathways/Results/HSVS_ER_CV", fold, 
		"_balanced.RData", sep=""))
}   

nsim=20000
start_nsim=10001
t0 <- proc.time()
## apply CV
snowfall::sfInit(parallel=TRUE, type="SOCK", cpus=5) 
snowfall::sfExport("X_s", "Y",  
  "g_index", "feature_types","types", "nsim", "start_nsim", "seed", 
  "HSVS_CV", "fold_id")
sfLibrary(msm)
snowfall::sfClusterApply(seq(1,5),function(x) 
  HSVS_CV(X_s, Y, g_index,
	  types, nsim, start_nsim, seed, fold_id, x))
snowfall::sfStop()
t1<-proc.time()-t0

t1 #70797.706 (19.7 hours)

