rm(list = ls())
setwd("/mnt/glusterfs/liz86/MOG_Regression/Exp_revision1/Simu1/Results")
#setwd("/mnt/glusterfs/liz86/MOG_Regression/Exp/simu1_v2")

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

summary_res("SOG_intercept_autoConv")


summary_res("Lasso_intercept")

summary_res("GEL")

load("Results/Veera.RData")
f_fdr_veera<-feature_FDR
f_for_veera<-feature_FOR
mse_veera<-MSE
auc_veera<-AUC

mean(mse_veera)
sd(mse_veera)/sqrt(length(f_fdr_sog))
mean(auc_veera)
sd(auc_veera)/sqrt(length(f_fdr_sog))

load("/mnt/glusterfs/liz86/MOG_Regression/Exp/simu1_v2/Results/Xu.RData")
f_fdr_xu<-feature_FDR
f_for_xu<-feature_FOR
mse_xu<-MSE
auc_xu<-AUC

mean(f_for_xu)
sd(f_for_xu)/sqrt(length(f_fdr_xu))
mean(f_fdr_xu)
sd(f_fdr_xu)/sqrt(length(f_fdr_xu))
mean(mse_xu)
sd(mse_xu)/sqrt(length(f_fdr_xu))
mean(auc_xu)
sd(auc_xu)/sqrt(length(f_fdr_xu))

load("/mnt/glusterfs/liz86/MOG_Regression/Exp/simu1_v2/Results/Lasso.RData")
f_fdr_lasso<-feature_FDR
f_for_lasso<-feature_FOR
mse_lasso<-MSE
auc_lasso<-AUC

mean(mse_lasso)
sd(mse_lasso)/sqrt(length(auc_lasso))
mean(auc_lasso)
sd(auc_lasso)/sqrt(length(auc_lasso))

load("Results/GroupLasso.RData")
f_fdr_gl<-feature_FDR
f_for_gl<-feature_FOR
mse_gl<-MSE
auc_gl<-AUC

mean(mse_gl)
sd(mse_gl)/sqrt(length(auc_lasso))
mean(auc_gl)
sd(auc_gl)/sqrt(length(auc_lasso))

load("Results/SGL.RData")
f_fdr_sgl<-feature_FDR
f_for_sgl<-feature_FOR
mse_sgl<-MSE
auc_sgl<-AUC

mean(mse_sgl)
sd(mse_sgl)/sqrt(length(auc_sgl))
mean(auc_sgl)
sd(auc_sgl)/sqrt(length(auc_sgl))



