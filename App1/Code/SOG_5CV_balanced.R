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

feature_types <- paste(featureNames_nodup, type_nodup, sep="_")

####### scale X ########
d_nodup <- d_nodup[which(ER %in% c("Positive", "Negative")), ]
ER <- ER[which(ER %in% c("Positive", "Negative"))]

X<-d_nodup
Y <- ER
Y <- (Y == "Positive")*1
F <- ncol(X)
feature_types <- sapply(1:ncol(X), function(i) paste(colnames(d_nodup)[i], 
    type_nodup[i], sep="_"))
X_s <- scale(X)

##### group index
g_index <- as.numeric(factor(featureNames_nodup))
types <- as.numeric(as.factor(type_nodup))

seed <- 20161120
save(X_s, Y, g_index, types, fold_id, feature_types, 
	burnInIter, totalIter, seed, file="Data/SOG_ER_CV_Pre_balanced_clean.RData")
save.image(file="/mnt/glusterfs/liz86/MOG_Regression/Exp/app1_8pathways/Data/SOG_ER_CV_Pre_balanced.RData")

## function for cross-validation
SOG_CV <- function(X_s, Y, g_index, 
	  types, totalIter, burnInIter, seed, fold_id, fold) {
	library(glmnet)
	sourceCpp("/mnt/glusterfs/liz86/MOG_Regression/Code/MOG_pkg/test_function/SOG_binary.cpp")
	source("/mnt/glusterfs/liz86/MOG_Regression/Code/MOG_pkg/test_function/SOG_binary.R") 
	X_s_train <- X_s[which(fold_id != fold), ]
	Y_train <- Y[which(fold_id != fold)]
	res <- SOG_binary(X_s_train, Y_train, 
		g_index, types=types, totalIter=totalIter, 
		burnInIter=burnInIter, seed=seed, lassoInit=FALSE)
	save(res, file=paste("/mnt/glusterfs/liz86/MOG_Regression/Exp/app1_8pathways/Results/SOG_ER_CV", fold, 
		"_balanced.RData", sep=""))
}

## apply CV
t0 <- proc.time()
snowfall::sfInit(parallel=TRUE, type="SOCK", cpus=5) 
snowfall::sfExport("X_s", "Y",  
  "g_index", "types", "totalIter", "burnInIter", "seed", 
  "SOG_CV", "fold_id")
sfLibrary(Rcpp)
sfLibrary(RcppEigen)
snowfall::sfClusterApply(seq(1,5),function(x) 
  SOG_CV(X_s, Y, g_index,
	  types, totalIter, burnInIter, seed, fold_id, x))
snowfall::sfStop()
t1<-proc.time()-t0

t1 #265.154 

