rm(list=ls())

setwd("/home06/liz86/BayesGL/exp/paper_simulation_3_v3_train_test")
se<-function(x)sd(x, na.rm=TRUE)/sqrt(length(x[!is.na(x)]))


load("Data/GroupLasso.RData")
F_FDR_GroupLasso<-F_FDR
F_FOR_GroupLasso<-F_FOR
MSE_GroupLasso<-MSE
AUC_GroupLasso<-AUC

apply(F_FDR,1,mean)
apply(F_FDR,1,se)
apply(F_FOR,1,mean)
apply(F_FOR,1,se)
apply(MSE,1,mean)
apply(MSE,1,se)
apply(AUC,1,mean)
apply(AUC,1,se)


load("Data/SparseGL.RData")
F_FDR_SGL<-F_FDR
F_FOR_SGL<-F_FOR
MSE_SGL<-MSE
AUC_SGL<-AUC
apply(MSE,1,mean)
apply(MSE,1,se)
apply(AUC,1,mean)
apply(AUC,1,se)


