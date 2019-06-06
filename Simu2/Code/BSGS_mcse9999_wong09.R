rm(list = ls())
library(pROC)
library(BSGS)
library(snowfall)
setwd("/mnt/glusterfs/liz86/MOG_Regression/Exp_revision1/Simu2/")
source("Code/BSGS_func.R")

simu2_BSGS(selected_tau2=1, mcse=9999, totalSimu=100, num_cpu=20, 
	burnInIter=20000, totalIter=30000)

simu2_BSGS(selected_tau2=5, mcse=9999, totalSimu=100, num_cpu=20, 
	burnInIter=20000, totalIter=30000)







