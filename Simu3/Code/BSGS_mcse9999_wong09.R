#########################################################
# Simulation 3 in paper (one-layer group structure)
# 
# 1/10/2018
# Li Zhu
########################################################
rm(list = ls())
library(pROC)
library(BSGS)
library(snowfall)
setwd("/mnt/glusterfs/liz86/MOG_Regression/Exp/simu3/")
source("Code/BSGS_func.R")

#simu3_BSGS(U=0.2, du=1, tau2_choice=seq(1,5), mcse=9999, 
#	totalSimu=100, num_cpu=10, 
#	burnInIter=100000, totalIter=200000)

simu3_BSGS(U=0.5, du=2, tau2_choice=seq(1,5), mcse=9999, 
	totalSimu=100, num_cpu=10, 
	burnInIter=100000, totalIter=200000)









