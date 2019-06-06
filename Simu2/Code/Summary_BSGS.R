#########################################################
# Simulation 1 in paper (one-layer group structure)
# 
# 1/6/2018
# Li Zhu
########################################################
rm(list = ls())
library(pROC)

setwd("/mnt/glusterfs/liz86/MOG_Regression/Exp_revision1/Simu2/")

summary_func <- function(res, mcse){
	print(paste("mean(AUC)=", round(mean(res$AUC[which(res$MCSE < mcse)],na.rm=TRUE), digits=2), sep=""))
	print(paste("SE(AUC)=", round(sd(res$AUC[which(res$MCSE < mcse)],na.rm=TRUE)/sqrt(sum(!is.na(res$AUC))), digits=2), sep=""))

	print(paste("mean(MSE)=", round(mean(res$MSE[which(res$MCSE < mcse)],na.rm=TRUE), digits=2), sep=""))
	print(paste("SE(MSE)=", round(sd(res$MSE[which(res$MCSE < mcse)],na.rm=TRUE)/sqrt(sum(!is.na(res$MSE))), digits=2), sep=""))

	print(paste("mean(FDR)=", round(mean(res$feature_FDR[which(res$MCSE < mcse)],na.rm=TRUE), digits=2),sep=""))
	print(paste("SE(FDR)=", round(sd(res$feature_FDR[which(res$MCSE < mcse)],na.rm=TRUE)/sqrt(sum(!is.na(res$feature_FDR))),digits=2), sep=""))

	print(paste("mean(FOR)=", round(mean(res$feature_FOR[which(res$MCSE < mcse)],na.rm=TRUE),digits=2), sep=""))
	print(paste("SE(FOR)=", round(sd(res$feature_FOR[which(res$MCSE < mcse)],na.rm=TRUE)/sqrt(sum(!is.na(res$feature_FOR))),digits=2), sep=""))
}

load("Results/BSGS_tau2_1_mcse_9999_res.RData")
res$delta_t #17802.073 5 hours (5 simulations per core, with 30k iterations)
sum(res$MCSE < 0.1) #10
summary_func(res, 0.1)
#[1] "mean(AUC)=0.87"
#[1] "SE(AUC)=0.02"
#[1] "mean(MSE)=23.66"
#[1] "SE(MSE)=3.1"
#[1] "mean(FDR)=0.26"
#[1] "SE(FDR)=0.06"
#[1] "mean(FOR)=0.07"
#[1] "SE(FOR)=0"

load("Results/BSGS_tau2_5_mcse_9999_res.RData")
res$delta_t #17802.073 
sum(res$MCSE < 0.1) #52
summary_func(res, 0.1)
#[1] "mean(AUC)=0.97"
#[1] "SE(AUC)=0"
#[1] "mean(MSE)=6.16"
#[1] "SE(MSE)=0.7"
#[1] "mean(FDR)=0.05"
#[1] "SE(FDR)=0.01"
#[1] "mean(FOR)=0.06"
#[1] "SE(FOR)=0"



