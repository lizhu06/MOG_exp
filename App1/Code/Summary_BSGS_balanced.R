###########################################
# Summary ER 5 cross-validation results
# Li Zhu
# 11/17/2016
###########################################
rm(list=ls())
library(sparcl)
library(gplots)
library(AUC)

setwd("/mnt/glusterfs/liz86/MOG_Regression/Exp/app1_8pathways/")

load("Data/SOG_ER_CV_Pre_balanced.RData")
nsim <- 20000
start_nsim <- 10001
res_cv <-list()
for(i in 1:5){
    load(paste("Results/BSGS_ER_CV", i, "_balanced.RData", sep=""))
    res_cv[[i]] <- res
}

AUC1_cv <- AUC2_cv <- rep(NA, 5)
for(i in 1:5){
    X_s_test <- X_s[fold_id == i, ]
    Y_test <- Y[fold_id == i]
    PROB <- matrix(NA, length(Y_test), ncol(res_cv[[1]]$coef))
    for(s in 1:ncol(res_cv[[i]]$coef)) {
        PROB[, s] <- c(pnorm(X_s_test[,order(g_index, decreasing=FALSE)] %*% 
            res_cv[[i]]$coef[,s] + 
            res_cv[[i]]$beta0_store[s]))
    }
    prob<-apply(PROB,1,mean)
    AUC2_cv[i] <- auc(roc(prob, factor(Y_test))) 

    ### plug-in posterior median of beta
    beta_median <- apply(res_cv[[i]]$coef, 1, median)
    prob<-pnorm(X_s_test[,order(g_index, decreasing=FALSE)]%*%beta_median + 
        median(res_cv[[i]]$beta0_store))
    AUC1_cv[i] <- auc(roc(prob, factor(Y_test)))  
}
mean(AUC1_cv) #0.940 (intercept: 0.941, balanced:0.947)
sd(AUC1_cv)/sqrt(5) #0.014 (intercept: 0.015, balanced: 0.009)

mean(AUC2_cv)  #0.944 (intercept: 0.945, balanced: 0.948)
sd(AUC2_cv)/sqrt(5) #0.011 (intercept: 0.014, balanced: 0.013)



################################
##### pathway impact score #####
################################
X_s_order <- X_s[, order(g_index, decreasing=FALSE)]
pm_per_pathway <- matrix(NA, C, 5)
FIS <- matrix(NA, length(unique(feature_types)), 5)
uni_features <- unique(feature_types[order(g_index, decreasing=FALSE)])
rownames(FIS) <- uni_features

for(i in 1:5){
    PM <- res_cv[[i]]$coef != 0

    ## sort features
    FIS[,i] <- pm <- apply(PM, 1, mean)

    ## sort pathways  
    pm_per_pathway[, i] <- sapply(1:C, function(x) mean(pm[which(colnames(X_s_order) %in% pathway_data[[x]])]))
    pathway_length <- sapply(1:C, function(x) sum(colnames(X_s_order) %in% pathway_data[[x]]))
    pm_per_feature <- sapply(1:C, function(x)
        paste(round(pm[which(colnames(X_s_order) %in% pathway_data[[x]])][order(pm[which(colnames(X_s_order) %in% pathway_data[[x]])], 
            decreasing=T)], 2), collapse="//"))
    feature_names <- sapply(1:C, function(x)
        paste(paste(colnames(X_s_order)[which(colnames(X_s_order) %in% 
            pathway_data[[x]])][order(pm[which(colnames(X_s_order) %in% 
                pathway_data[[x]])], decreasing=T)], 
            type_nodup[which(colnames(X_s_order) %in% pathway_data[[x]])][order(pm[which(colnames(X_s_order) %in% pathway_data[[x]])], 
            decreasing=T)], sep="_"), collapse="//"))

    pathway_impact <- cbind(names(pathway_data), round(pm_per_pathway[,i], 3),
      pathway_length, pm_per_feature, feature_names)[order(pm_per_pathway[,i], decreasing=T), ]
    colnames(pathway_impact) <- c("pathway_name", "impact_score", 
        "pathway_size", "PM_features", "feature_names")
    write.csv(pathway_impact, 
        file=paste("Results/BSGS_ER_pathway_impact", i ,
            "_balanced.csv", sep=""), row.names=FALSE)
}

## sort pathways
rownames(pm_per_pathway) <- names(pathway_data)
mean_PIS <- apply(pm_per_pathway, 1, mean)
names(mean_PIS) <- names(pathway_data)
mean_PIS_sort <- mean_PIS[order(mean_PIS, decreasing=TRUE)]
mean_PIS_sort[1:5]
#Estrogen signaling
# 0.020

## sort features
mean_FIS <- apply(FIS, 1, mean)
names(mean_FIS) <- rownames(FIS)
mean_FIS_sort <- mean_FIS[order(mean_FIS, decreasing=TRUE)]
mean_FIS_sort[1:5]
#  ESR1_gene   NME3_gene  ADCY9_gene ESR1_methyl   SHC2_gene 
#    1.00000     0.43764     0.14978     0.14590     0.12398 

## sort dup features
index <- which(names(pathway_data)=="Estrogen signaling pathway_Homo sapiens_hsa04915")
mean_dfis <- apply(FIS, 1, mean)
uni_mean_dfis <- unique(mean_dfis)
uni_mean_dfis_sort <- uni_mean_dfis[order(uni_mean_dfis, decreasing=TRUE)]
num_nonzero_BSGS <- sapply(1:length(uni_mean_dfis_sort), function(x) 
    sum(mean_dfis >= uni_mean_dfis_sort[x]))
num_nonzero_ER_BSGS <- sapply(1:length(uni_mean_dfis_sort), function(x) 
    sum(colnames(X_s_order)[mean_dfis >= uni_mean_dfis_sort[x]] %in% pathway_data[[index]]))

save(num_nonzero_BSGS, file="Results/num_nonzero_BSGS_balanced.RData")
save(num_nonzero_ER_BSGS, file="Results/num_nonzero_ER_BSGS_balanced.RData")




