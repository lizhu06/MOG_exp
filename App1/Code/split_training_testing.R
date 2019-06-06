rm(list=ls())

##### duplicate genes belonging to multiple groups
setwd("/mnt/glusterfs/liz86/MOG_Regression/Exp_revision1/App1/")
#load("ER.RData")
load("Data/d_nodup.RData")
load("Data/ER.RData")
load("Data/type_nodup.RData")
load("Data/pathway_data_no_singleton.RData")

all(names(ER) == rownames(d_nodup))
X <- d_nodup
X_s <- scale(X)
Y <- ER
Y <- Y=="Positive"
feature_types <- sapply(1:ncol(X), function(i) paste( 
	colnames(d_nodup)[i],type_nodup[i],sep="_"))
types <- as.numeric(as.factor(type_nodup))

## create U1 and U2
uni_genes <- unique(colnames(d_nodup))
U1 <- matrix(0, ncol(X_s), length(uni_genes))
colnames(U1) <- uni_genes
for(i in 1:ncol(X_s)){
	U1[i, colnames(X_s)[i]] <- 1
}
U2 <- matrix(0, length(uni_genes), length(pathway_data))
rownames(U2) <- uni_genes
for(i in 1:length(pathway_data)){
	U2[pathway_data[[i]], i]<-1
}

### split data for training and testing
set.seed(20161120)
fold_id <- rep(NA, length(Y))
sum(Y==1) #560
sum(Y==0) #167
fold_id[Y==1] <- sample(c(rep(1, 112), rep(2, 112), rep(3, 112), rep(4, 112), rep(5, 112)))
fold_id[Y==0] <- sample(c(rep(1, 34), rep(2, 34), rep(3, 33), rep(4, 33), rep(5, 33)))
for(i in 1:5){
	print(mean(Y[fold_id==i]))
}

save(X_s, Y, U1, U2, types, fold_id, feature_types, 
	file="Data/Data_clean.RData")
