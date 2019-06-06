########################
# BSGS-SS predict ER
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


## function for cross-validation
BSGS_CV <- function(X_s, Y, g_index, 
	  feature_types, types, nsim, start_nsim, seed, fold_id, fold) {
	source("/mnt/glusterfs/liz86/MOG_Regression/R/BSGSSS/BSGSSS_pkg_binary_v2.R")  
	source("/mnt/glusterfs/liz86/MOG_Regression/R/BSGSSS/MBSGS-internal_binary_v2.R")  
	X_s_train <- X_s[which(fold_id != fold), ]
	Y_train <- Y[which(fold_id != fold)]
	X_s_train_sort <- X_s_train[, order(g_index, decreasing=FALSE)]
	feature_types_sort <- feature_types[order(g_index, decreasing=FALSE)]
	gsize <- sapply(1:max(g_index), function(x) sum(g_index==x))
	res <- BSGSSS_binary(Y_train, X_s_train_sort, niter=nsim, burnin=start_nsim-1, group_size=gsize, num_update = 100, niter.update = 100)  # niter is the total iterations, final iteration is niter-burnin
	save(res, file=paste("/mnt/glusterfs/liz86/MOG_Regression/Exp/app1_8pathways/Results/BSGS_ER_CV", fold, 
		"_balanced.RData", sep=""))
}

t0 <- proc.time()
## apply CV

nsim <- 20000
start_nsim <- 10001
snowfall::sfInit(parallel=TRUE, type="SOCK", cpus=5) 
snowfall::sfExport("X_s", "Y",  
  "g_index", "feature_types","types", "nsim", "start_nsim", "seed", 
  "BSGS_CV", "fold_id")
snowfall::sfClusterApply(seq(1,5),function(x) 
  BSGS_CV(X_s, Y, g_index,
	  feature_types, types, nsim, start_nsim, seed, fold_id, x))
snowfall::sfStop()
t1<-proc.time()-t0

t1 #4910.027 (1.4 hours)
