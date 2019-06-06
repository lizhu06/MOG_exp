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
    load(paste("Results/HSVS_ER_CV", i, "_balanced.RData", sep=""))
    res_cv[[i]] <- res
}

AUC1_cv <- AUC2_cv <- rep(NA, 5)
for(i in 1:5){
    X_s_test <- X_s[fold_id == i, ]
    Y_test <- Y[fold_id == i]
    PROB <- matrix(NA, length(Y_test), (nsim - start_nsim) + 1)
    for(s in 1:(nsim - start_nsim + 1)) {
        PROB[, s] <- c(pnorm(X_s_test %*% res_cv[[i]][[1]][,s]))
    }
    prob<-apply(PROB,1,mean)
    AUC2_cv[i] <- auc(roc(prob, factor(Y_test))) 

    ### plug-in posterior median of beta
    beta_median <- apply(res_cv[[i]][[1]], 1, median)
    prob<-pnorm(X_s_test%*%beta_median)
    AUC1_cv[i] <- auc(roc(prob, factor(Y_test)))  
}
mean(AUC1_cv) #0.942 (balanced: 0.942)
sd(AUC1_cv)/sqrt(5) #0.014 (balanced: 0.012)

mean(AUC2_cv)  #0.945 (balanced: 0.944)
sd(AUC2_cv)/sqrt(5) #0.013 (balanced: 0.012)



################################
##### pathway impact score #####
################################
pm_per_pathway <- matrix(NA, C, 5)
FIS <- matrix(NA, length(unique(feature_types)), 5)
uni_features <- unique(feature_types)
rownames(FIS) <- uni_features

for(i in 1:5){
    PM <- res_cv[[i]][[1]] != 0

    ## sort features
    FIS[,i] <- pm <- apply(PM, 1, mean)

    ## sort pathways  
    pm_per_pathway[, i] <- sapply(1:C, function(x) mean(pm[which(colnames(X_s) %in% pathway_data[[x]])]))
    pathway_length <- sapply(1:C, function(x) sum(colnames(X_s) %in% pathway_data[[x]]))
    pm_per_feature <- sapply(1:C, function(x)
        paste(round(pm[which(colnames(X_s) %in% pathway_data[[x]])][order(pm[which(colnames(X_s) %in% pathway_data[[x]])], decreasing=T)], 
            2), collapse="//"))
    feature_names <- sapply(1:C, function(x)
        paste(paste(colnames(X_s)[which(colnames(X_s) %in% pathway_data[[x]])][
            order(pm[which(colnames(X_s) %in% pathway_data[[x]])], 
                decreasing=T)], 
            type_nodup[which(colnames(X_s) %in% pathway_data[[x]])][order(pm[which(colnames(X_s) %in% pathway_data[[x]])], 
            decreasing=T)], sep="_"), collapse="//"))

    pathway_impact <- cbind(names(pathway_data), round(pm_per_pathway[,i], 3),
      pathway_length, pm_per_feature, feature_names)[order(pm_per_pathway[,i], decreasing=T), ]
    colnames(pathway_impact) <- c("pathway_name", "impact_score", 
        "pathway_size", "PM_features", "feature_names")
    write.csv(pathway_impact, 
        file=paste("Results/HSVS_ER_pathway_CV", i ,
            "_balanced.csv", sep=""), row.names=FALSE)
}

## sort pathways
rownames(pm_per_pathway) <- names(pathway_data)
mean_PIS <- apply(pm_per_pathway, 1, mean)
names(mean_PIS) <- names(pathway_data)
mean_PIS_sort <- mean_PIS[order(mean_PIS, decreasing=TRUE)]
mean_PIS_sort[1:5]
#Estrogen signaling
#0.027

## sort features
mean_FIS <- apply(FIS, 1, mean)
names(mean_FIS) <- rownames(FIS)
mean_FIS_sort <- mean_FIS[order(mean_FIS, decreasing=TRUE)]
mean_FIS_sort[1:5]
#  ESR1_gene    ESR1_CNV ESR1_methyl   NME3_gene    NME3_CNV 
#     1.0000      1.0000      1.0000      0.3084      0.3084 

## sort dup features
index <- which(names(pathway_data)=="Estrogen signaling pathway_Homo sapiens_hsa04915")
mean_dfis <- apply(FIS, 1, mean)
uni_mean_dfis <- unique(mean_dfis)
uni_mean_dfis_sort <- uni_mean_dfis[order(uni_mean_dfis, decreasing=TRUE)]
num_nonzero_HSVS <- sapply(1:length(uni_mean_dfis_sort), function(x) 
    sum(mean_dfis >= uni_mean_dfis_sort[x]))
num_nonzero_ER_HSVS <- sapply(1:length(uni_mean_dfis_sort), function(x) 
    sum(colnames(X_s)[mean_dfis >= uni_mean_dfis_sort[x]] %in% pathway_data[[index]]))

save(num_nonzero_HSVS, file="Results/num_nonzero_HSVS_balanced.RData")
save(num_nonzero_ER_HSVS, file="Results/num_nonzero_ER_HSVS_balanced.RData")

