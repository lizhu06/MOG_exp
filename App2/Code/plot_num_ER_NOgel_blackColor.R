rm(list=ls())
##### duplicate genes belonging to multiple groups
setwd("/mnt/glusterfs/liz86/MOG_Regression/Exp_revision1/App2/")

old_dir <- "/mnt/glusterfs/liz86/MOG_Regression/Exp/app1_morePathways/"

### plot all variable selection
load("Results/num_nonzero_MOG_balanced_nodup_lassoInit.RData")
load("Results/num_nonzero_ER_MOG_balanced_nodup_lassoInit.RData")
nonzero_MOG <- num_nonzero_MOG_nodup
nonzero_ER_MOG <- num_nonzero_ER_MOG_nodup

load(paste0(old_dir, "Results/num_nonzero_SOG_balanced.RData"))
load(paste0(old_dir, "Results/num_nonzero_ER_SOG_balanced.RData"))
nonzero_SOG <- num_nonzero_SOG
nonzero_ER_SOG <- num_nonzero_ER_SOG

for(i in 1:5){
	load(paste0(old_dir, "Results/Lasso_nonzero_balanced.RData"))
	load(paste0(old_dir, "Results/Lasso_nonzero_ER_balanced.RData"))
	nonzero_lasso <- nonzero[i,]
	nonzero_ER_lasso <- nonzero_ER[i,]

	load(paste0(old_dir, "Lasso_nonzero_GL_balanced.RData"))
	load(paste0(old_dir, "Lasso_nonzero_ER_GL_balanced.RData"))
	nonzero_GL <- nonzero_GL[i,]
	nonzero_ER_GL <- nonzero_ER_GL[i,]

	load(paste0(old_dir, "nonzero_SGL_balanced.RData"))
	load(paste0(old_dir, "nonzero_ER_SGL_balanced.RData"))
	nonzero_SGL <- nonzero_SGL[i,]
	nonzero_ER_SGL <- nonzero_ER_SGL[i,]

	load("Results/nonzero_GEL_balanced.RData")
	load("Results/nonzero_ER_GEL_balanced.RData")
	nonzero_GEL <- nonzero_GEL[i,]
	nonzero_ER_GEL <- nonzero_ER_GEL[i,]

	load("Results/nonzero_TSG_balanced.RData")
	load("Results/nonzero_ER_TSG_balanced.RData")
	nonzero_TGL <- nonzero_TSG[i,]
	nonzero_ER_TGL <- nonzero_ER_TSG[i,]


	pdf(paste("Results/ER_genesInER_cv", i,"_balanced_noGEL_blackColor.pdf",sep=""),width=5,height=5)
	plot(nonzero_MOG,nonzero_ER_MOG,xlim=c(1,200),ylim=c(1,150),
		ylab="Number of features in ER pathway",xlab="Number of features selected",
		col=1,type="b",pch=1, cex=0.5)
	lines(nonzero_SOG,nonzero_ER_SOG,col=1,type="b",pch=2, cex=0.5)
	lines(nonzero_lasso,nonzero_ER_lasso,col=1,type="b",pch=3, cex=0.5)
	lines(nonzero_GL,nonzero_ER_GL,col=1,type="b",pch=4, cex=0.5)
	lines(nonzero_SGL,nonzero_ER_SGL,col=1,type="b",pch=5, cex=0.5)
	#lines(nonzero_GEL,nonzero_ER_GEL,col=6,type="b",pch=6)
	lines(nonzero_TGL,nonzero_ER_TGL,col=1,type="b",pch=6, cex=0.5)
	abline(0,1,lty=2)
	#abline(0,0.11,lty=2)
	legend("topleft",c("MOG","SOG","Lasso","GL","SGL",  "TGL"),col=1,
		lty=1,pch=1:6)
	dev.off()
}
