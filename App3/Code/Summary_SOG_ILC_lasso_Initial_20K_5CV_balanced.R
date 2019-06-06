###########################################
# Summary ER 5 cross-validation results
# Li Zhu
# 11/17/2016
###########################################
rm(list=ls())
library(sparcl)
library(gplots)
library(AUC)

setwd("/mnt/glusterfs/liz86/MOG_Regression/Exp/app2/")

load("Data/SOG_ILC_CV_Pre_balanced.RData")
res_cv <-list()
for(i in 1:5){
    load(paste("Results/SOG_ILC_CV", i, "_balanced.RData", sep=""))
    res_cv[[i]] <- res
}

AUC1_cv <- AUC2_cv <- rep(NA, 5)
for(i in 1:5){
    X_s_test <- X_s[fold_id == i, ]
    Y_test <- Y[fold_id == i]
    PROB <- matrix(NA, length(Y_test), totalIter-burnInIter)
    for(s in 1:(totalIter-burnInIter)) {
        PROB[, s] <- c(pnorm(X_s_test %*% res_cv[[i]]$BETA[-1,s+burnInIter]+
            res_cv[[i]]$BETA[1,s+burnInIter]))
    }
    prob<-apply(PROB,1,mean)
    AUC2_cv[i] <- auc(roc(prob, factor(Y_test))) 

    ### plug-in posterior median of beta
    beta_median <- apply(res_cv[[i]]$BETA[,-seq(1,burnInIter)], 1, median)
    prob<-pnorm(X_s_test%*%beta_median[-1]+beta_median[1])
    AUC1_cv[i] <- auc(roc(prob, factor(Y_test)))  
}
mean(AUC1_cv) #0.913 (0.922, balanced: 0.911)
sd(AUC1_cv)/sqrt(5) #0.013 (0.006, balanced: 0.010)

mean(AUC2_cv)  #0.955 (0.942, balanced: 0.950)
sd(AUC2_cv)/sqrt(5) #0.003 (0.007, balanced: 0.012)



################################
##### pathway impact score #####
################################
pm_per_pathway <- matrix(NA, C, 5)
FIS <- matrix(NA, length(unique(feature_types)), 5)
uni_features <- unique(feature_types)
rownames(FIS) <- uni_features

for(i in 1:5){
    PM <- res_cv[[i]]$BETA[-1, -seq(1,burnInIter)] != 0

    ## sort features
    FIS[,i] <- pm <- apply(PM, 1, mean)

    ## sort pathways  
    pm_per_pathway[, i] <- sapply(1:C, function(x) mean(pm[c_index == x]))
    pathway_length <- sapply(1:C, function(x) sum(c_index == x))
    pm_per_feature <- sapply(1:C, function(x)
        paste(round(pm[c_index == x][order(pm[c_index == x], decreasing=T)], 
            2), collapse="//"))
    feature_names <- sapply(1:C, function(x)
        paste(paste(uni_genes[true_gene_ind_dup[c_index == x]][
            order(pm[c_index == x], decreasing=T)], 
            type_dup[c_index == x][order(pm[c_index == x], 
            decreasing=T)], sep="_"), collapse="//"))

    pathway_impact <- cbind(names(pathway_data), round(pm_per_pathway[,i], 3),
      pathway_length, pm_per_feature, feature_names)[order(pm_per_pathway[,i], decreasing=T), ]
    colnames(pathway_impact) <- c("pathway_name", "impact_score", 
        "pathway_size", "PM_features", "feature_names")
    write.csv(pathway_impact, 
        file=paste("Results/SOG_ILC_pathway_impact_LassoInitial_20k_CV", i ,
            "_balanced.csv", sep=""), row.names=FALSE)
}

## sort pathways
rownames(pm_per_pathway) <- names(pathway_data)
mean_PIS <- apply(pm_per_pathway, 1, mean)
names(mean_PIS) <- names(pathway_data)
mean_PIS_sort <- mean_PIS[order(mean_PIS, decreasing=TRUE)]
mean_PIS_sort[1:5] #Viral myocarditis, 0.044

## sort features
mean_FIS <- apply(FIS, 1, mean)
names(mean_FIS) <- rownames(FIS)
mean_FIS_sort <- mean_FIS[order(mean_FIS, decreasing=TRUE)]
mean_FIS_sort[1:5]
#   CDH1_gene  MAP3K1_gene SHROOM1_gene   LAMA3_gene ALDH1B1_gene 
#     1.00000      0.62300      0.51688      0.51446      0.44064 
