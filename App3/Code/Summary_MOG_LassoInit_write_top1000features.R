rm(list=ls())
#library(sparcl)
#library(gplots)
library(pROC)

setwd("/mnt/glusterfs/liz86/MOG_Regression/Exp_revision1/App3/")
load("Data/pathway_data_filtered.RData")
load("Data/Data_clean.RData")  

res_cv <-list()
for(i in 1:5){
    load(paste("Results/MOG_ER_CV", i, "_balanced_LassoInit.RData", sep="")) # A mistake in file name, should be Histology not ER
    res_cv[[i]] <- res
}

FIS <- matrix(NA, length(unique(feature_types)), 5) # unique features
uni_features <- unique(feature_types)
rownames(FIS) <- uni_features

for(i in 1:5){
    # selection probability (already no-duplicated)
    PM <- res_cv[[i]]$BETA[-1, ] != 0
    FIS[,i] <- pm <- apply(PM, 1, mean)
    names(pm) <- colnames(X_s)
}

## sort features
mean_FIS <- apply(FIS, 1, mean)
names(mean_FIS) <- rownames(FIS)
mean_FIS_sort <- mean_FIS[order(mean_FIS, decreasing=TRUE)]
mean_FIS_sort[1:5]

num_top_gene <- 500
FIS_table <- matrix(NA, num_top_gene, 3)
for(i in 1:num_top_gene){
    gene_name <- strsplit(names(mean_FIS_sort)[i], split="_")[[1]][1]
    pathway_index <- sapply(1:length(pathway_data), function(x) gene_name %in% 
        pathway_data[[x]])
    in_pathway_name <- names(pathway_data)[pathway_index]
    FIS_table[i, 1] <- names(mean_FIS_sort)[i]
    FIS_table[i, 2] <- mean_FIS_sort[i]
    FIS_table[i, 3] <- paste(in_pathway_name, collapse="//")
    print(i)
}
colnames(FIS_table) <- c("FeatureName", "FIS", "Pathways")
write.csv(FIS_table, file="Results/FIS_table_top500.csv")
