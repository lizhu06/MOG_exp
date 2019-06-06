#########################################################
# Simulation 1 in paper (one-layer group structure)
# 
# 1/6/2018
# Li Zhu
########################################################
rm(list = ls())
library(pROC)

setwd("/mnt/glusterfs/liz86/MOG_Regression/Exp/simu3/")

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

load("Results/BSGS_U_0.2_mcse_9999_res.RData")
res$delta_t #39822.599 (11 hours)
sum(res$MCSE < 0.1) #97
summary_func(res, 0.1)

#[1] "mean(AUC)=0.86"
#[1] "SE(AUC)=0"
#[1] "mean(MSE)=29.64"
#[1] "SE(MSE)=1.05"
#[1] "mean(FDR)=0.03"
#[1] "SE(FDR)=0.01"
#[1] "mean(FOR)=0.25"
#[1] "SE(FOR)=0"

load("Results/BSGS_U_0.5_mcse_9999_res.RData")
res$delta_t  # 41854.726 (11.63 hours)
sum(res$MCSE < 0.1) #9
sum(res$MCSE < 0.5) # 77
summary_func(res, 0.1)
#[1] "mean(AUC)=0.94"
#[1] "SE(AUC)=0.01"
#[1] "mean(MSE)=17.01"
#[1] "SE(MSE)=1.36"
#[1] "mean(FDR)=0.08"
#[1] "SE(FDR)=0.02"
#[1] "mean(FOR)=0.04"
#[1] "SE(FOR)=0"




