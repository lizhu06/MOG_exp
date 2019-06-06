rm(list=ls())
library(glmnet)
library(pROC)
library(RcppEigen)
library(Rcpp)
library(coda)
library(snowfall)

##### duplicate genes belonging to multiple groups
setwd("/mnt/glusterfs/liz86/MOG_Regression/Exp_revision1/App1/")
load("Data/Data_clean.RData")  

## function for cross-validation
MOG_CV <- function(X_s, Y, U1, U2, 
	  types, burnInIter, keepIter, maxIter, seed, fold_id, fold,
	  lassoInit) {
	library(glmnet)
	library(coda)
	sourceCpp("/mnt/glusterfs/liz86/MOG_Regression/Code/MOG_revision1/src/MOG_binary_new.cpp")
	source("/mnt/glusterfs/liz86/MOG_Regression/Code/MOG_revision1/R/MOG_binary_new.R") 
	X_s_train <- X_s[which(fold_id != fold), ]
	Y_train <- Y[which(fold_id != fold)]
	
	res <- MOG_binary(X_s_train, Y_train, U1, U2, types=types, 
  lassoInit=lassoInit,  weight_Tj=TRUE, weight_Dk=TRUE, weight_Rj=FALSE,
  pi2_prop_n=10,  pi3_prop_n=10, MH_ind=1,
  burnInIter=burnInIter, keepIter=keepIter, maxIter=maxIter, 
  printInt=1000, seed=seed, 
  alpha1=1, beta1=1, alpha2=1, beta2=1, alpha3=1, beta3=1, 
  weight_s2=FALSE)
	
	save(res, file=paste("/mnt/glusterfs/liz86/MOG_Regression/Exp_revision1/App1/Results/MOG_ER_CV", fold, 
		"_balanced_LassoInit.RData", sep=""))
}

burnInIter <- 10000
keepIter <- 10000
maxIter <- 50000
seed <- 20161120
lassoInit <- TRUE
## apply CV
t0 <- proc.time()
snowfall::sfInit(parallel=TRUE, type="SOCK", cpus=5) 
snowfall::sfExport("X_s", "Y",  "U1", "U2",
  "types", "burnInIter", "keepIter", "maxIter", "seed", 
  "MOG_CV", "fold_id", "lassoInit")
sfLibrary(Rcpp)
sfLibrary(RcppEigen)
snowfall::sfClusterApply(seq(1,5),function(x) 
  MOG_CV(X_s, Y, U1, U2, 
	  types, burnInIter, keepIter, maxIter, seed, fold_id, x, lassoInit))
snowfall::sfStop()
t1<-proc.time()-t0



