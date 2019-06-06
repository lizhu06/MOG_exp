rm(list=ls())

##### duplicate genes belonging to multiple groups
setwd("/mnt/glusterfs/liz86/MOG_Regression/Exp_revision1/App3/")
#load("ER.RData")
load("Data/d_nodup.RData")
load("Data/Hist.RData")
load("Data/type_nodup.RData")
load("Data/pathway_data_filtered.RData")

all(names(Hist) == rownames(d_nodup))
X <- d_nodup
X_s <- scale(X)
Y <- Hist
Y <- Y=="ILC"
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
set.seed(20161116)
fold_id <- rep(NA, length(Y))
sum(Y==1) #173
sum(Y==0) #496
fold_id[Y==1] <- sample(c(rep(1, 35), rep(2, 35), rep(3, 35), rep(4, 34), rep(5, 34)))
fold_id[Y==0] <- sample(c(rep(1, 100), rep(2, 99), rep(3, 99), rep(4, 99), rep(5, 99)))
for(i in 1:5){
	print(mean(Y[fold_id==i]))
}

save(X_s, Y, U1, U2, types, fold_id, feature_types, 
	file="Data/Data_clean.RData")
