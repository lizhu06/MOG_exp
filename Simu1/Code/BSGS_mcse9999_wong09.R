#########################################################
# Simulation 1 in paper (one-layer group structure)
# 
# 1/6/2018
# Li Zhu
########################################################
rm(list = ls())
library(pROC)
library(BSGS)
library(snowfall)
setwd("/mnt/glusterfs/liz86/MOG_Regression/Exp/simu1_v2/")
source("Code/BSGS_func.R")

simu1_BSGS(selected_tau2=1, mcse=9999, num_cpu=20, totalSimu=100, 
	burnInIter=20000, totalIter=30000)

simu1_BSGS(selected_tau2=5, mcse=9999, num_cpu=20, totalSimu=100, 
	burnInIter=20000, totalIter=30000)







