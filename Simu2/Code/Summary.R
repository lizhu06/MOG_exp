rm(list = ls())
setwd("/mnt/glusterfs/liz86/MOG_Regression/Exp_revision1/Simu2/Results")

summary_res <- function(res_name){
	file_name <- paste0(res_name, ".RData")
	load(file_name)

	f_fdr_sog <- feature_FDR
	f_for_sog<-feature_FOR
	mse_sog<-MSE
	auc_sog<-AUC

	cat(paste0("AUC mean=", round(mean(auc_sog, na.rm=TRUE), digits=2), "\n"))
	cat(paste0("AUC SE=", round(sd(auc_sog, na.rm=TRUE)/
		sqrt(sum(!(is.na(auc_sog))))
		, digits=2), "\n"))

	cat(paste0("MSE mean=", round(mean(mse_sog, na.rm=TRUE), digits=2), "\n"))
	cat(paste0("MSE SE=", round(sd(mse_sog, na.rm=TRUE)/
		sqrt(sum(!(is.na(mse_sog))))
		, digits=3), "\n"))

	cat(paste0("FDR mean= ", round(mean(f_fdr_sog, na.rm=TRUE), digits=2), "\n"))
	cat(paste0("FDR SE= ", round(sd(f_fdr_sog, na.rm=TRUE)/
		sqrt(sum(!is.na(f_fdr_sog)))
		, digits=2), "\n"))
	

	cat(paste0("FOR mean=", round(mean(f_for_sog, na.rm=TRUE), digits=2), "\n"))
	cat(paste0("FOR SE=", round(sd(f_for_sog, na.rm=TRUE)/
		sqrt(sum(!is.na(f_for_sog)))
		, digits=2), "\n"))
	
}

summary_res("SOG_MH_intercept_autoConv")
summary_res("GEL")


#load("Xu.RData")
summary_res("Xu")

#load("Veera.RData")
summary_res("Veera")