########################################
# preprocess data for predict ER+/-
# 	using only 8 KEGG pathways
# Li Zhu
# 12/09/2016
#########################################
rm(list=ls())
setwd("/home06/liz86/BayesGL/exp/paper_app1_8pathways/Data")
load("/home06/liz86/BayesGL/Data/TCGA/expr_meth_cnv.RData")
load("/home06/liz86/BayesGL/Data/TCGA/omicsTyps.RData")
load("/home06/liz86/BayesGL/Data/TCGA/featureNames.RData")

### transform methylation to M-value
d2 <- d
# transfrom to log(odds) for methylation
for(n in 1:nrow(d)){
	d2[n,which(omicsTyps=="methyl")]<-log(d[n,which(omicsTyps=="methyl")]/
		(1-d[n,which(omicsTyps=="methyl")]))
}

### patients ####
patients<-rownames(d2)
pat_names<-substr(patients,start=1,stop=12)
rownames(d2) <- pat_names
load("/home05/liz86/Steffi/Kevin_IDC_ILC_DE/Data/tcgaClinical_10_14_15_Processed.Rda")
clin<-tcga_clinical_10_14_15
sum(pat_names %in% rownames(clin))  #770
ER<-as.character(clin[match(pat_names,rownames(clin)),"er_status_by_ihc"])
names(ER) <- pat_names

# exclude outcome
d3<-d2[which(ER %in%c("Positive","Negative")),]
ER<-ER[which(ER %in%c("Positive","Negative"))]
all(names(ER) == rownames(d3))

featureNames3<-featureNames
omicsTyps3<-omicsTyps

dim(d3)  #727 14976
length(featureNames3)  #14976
length(omicsTyps3)   #14976


if(F){
	## pathway dataset
library('GSA')
## pathway dataset from Joey
load("/home06/liz86/MetDiffNet/package/MetaDiffNetwork/data/pathway.database.v4.0.Rdata")
num_gene_path<-sapply(1:length(pathway.database),function(x)length(pathway.database[[x]]))
path_name<-sapply(1:length(names(pathway.database)),function(x)strsplit(names(pathway.database)[x],split="_")[[1]][1])
initial_name<-sapply(1:length(names(pathway.database)),function(x)strsplit(names(pathway.database)[x],split="_")[[1]][2])

names(pathway.database)[which(path_name=="GO" & initial_name=="POSITIVE")]

length(pathway.database[num_gene_path2<=50 & path_name=="REACTOME"])  #478
length(pathway.database[num_gene_path2<=50 & path_name=="KEGG"])  #87
length(pathway.database[num_gene_path2<=50 & path_name=="BIOCARTA"])  #213
length(pathway.database[num_gene_path2<=50 & path_name=="GO"]) #994

length(pathway.database[path_name=="REACTOME"])  #674
length(pathway.database[path_name=="KEGG"])      #186
length(pathway.database[path_name=="BIOCARTA"])  #217
length(pathway.database[path_name=="GO"])        #1454

}

## load all KEGG pathway
load("/home06/liz86/BayesGL/Data/PathwayDB/KEGG_2016.RData")
num_gene_path<-sapply(1:length(KEGG_2016),function(x)length(KEGG_2016[[x]]))
path_name<-sapply(1:length(names(KEGG_2016)),function(x)strsplit(names(KEGG_2016)[x],split="_")[[1]][1])
"Estrogen signaling pathway" %in% path_name

database3<-KEGG_2016
length(database3)  #293

uni_gene<-unique(featureNames3) # unique genes in dataset
## exclude genes not in dataset
database4<-lapply(1:length(database3),function(x)intersect(database3[[x]],uni_gene))
names(database4)<-names(database3)

## Calculate gene numbers in each pathway
pathway_length<-sapply(1:length(database4),function(x)length(database4[[x]]))
pathway_length
sum(pathway_length<=50 & pathway_length>=40)  #26

## t-test for genes
pathway_gene<-unique(unlist(database4)) # unique genes in both dataset and pathway database
length(pathway_gene)   #1978
ttest_p<-sapply(1:length(pathway_gene),function(x)
	t.test(d3[which(ER=="Positive"),pathway_gene[x]],
		d3[which(ER=="Negative"),pathway_gene[x]])$p.value)
ttest_padj<-p.adjust(ttest_p,method="BH")
names(ttest_padj)<-pathway_gene
prop_sig_genes<-sapply(1:length(database4),function(x) 
	sum(ttest_padj[database4[[x]]]<0.05)/(length(database4[[x]])+1))

## select a subset of pathways
pathway_data<-database4[which(pathway_length<=50 & 
	pathway_length>=40 & prop_sig_genes>=0.8)]
length(pathway_data)  #8
uni_genes<-unique(unlist(pathway_data))
length(uni_genes)  #276
sum(featureNames3 %in% uni_genes) #824
path_length<-sapply(1:length(pathway_data),
	function(x)length(pathway_data[[x]]))

## data without duplication, only keep pathway genes
d_nodup<-d3[,which(featureNames3 %in% uni_genes)]
featureNames_nodup<-featureNames3[which(featureNames3 %in% uni_genes)]
type_nodup<-omicsTyps3[which(featureNames3 %in% uni_genes)]
dim(d_nodup) #727 824
length(featureNames_nodup) #824
length(type_nodup) #824


save(d_nodup,file="d_nodup.RData")
save(featureNames_nodup,file="featureNames_nodup.RData")
save(type_nodup,file="type_nodup.RData")

## duplication
d_dup<-rep(NA,nrow(d3))
pathway_ind_dup<-NA
fake_gene_ind_dup<-NA
true_gene_ind_dup<-NA
g_ind<-0
type_dup<-NA
for(i in 1:length(pathway_data)){
	for(j in 1:length(pathway_data[[i]])){
		g_ind<-g_ind+1
		d_dup<-cbind(d_dup,d3[,which(featureNames3==pathway_data[[i]][j])])
		fake_gene_ind_dup<-c(fake_gene_ind_dup,rep(g_ind,sum(featureNames3==pathway_data[[i]][j])))
		true_gene_ind_dup<-c(true_gene_ind_dup,rep(which(uni_genes==pathway_data[[i]][j]),sum(featureNames3==pathway_data[[i]][j])))
		pathway_ind_dup<-c(pathway_ind_dup,rep(i,sum(featureNames3==pathway_data[[i]][j])))
		type_dup<-c(type_dup,omicsTyps3[featureNames3==pathway_data[[i]][j]])
	}
	print(i)
}
d_dup<-d_dup[,-1]
pathway_ind_dup<-pathway_ind_dup[-1]
fake_gene_ind_dup<-fake_gene_ind_dup[-1]
true_gene_ind_dup<-true_gene_ind_dup[-1]
type_dup<-type_dup[-1]

max(pathway_ind_dup)  #8
max(true_gene_ind_dup) #276
max(fake_gene_ind_dup) #371

dim(d_dup)#727 1109
length(pathway_ind_dup)    #1109
length(fake_gene_ind_dup)  #1109
length(true_gene_ind_dup)  #1109
length(type_dup) #1109

## Save Data
save(d_dup,file="d_dup_no_singleton.RData")
save(uni_genes,file="uni_genes_no_singleton.RData")
save(pathway_ind_dup,file="pathway_ind_dup_no_singleton.RData")
save(fake_gene_ind_dup,file="fake_gene_ind_dup_no_singleton.RData")
save(true_gene_ind_dup,file="true_gene_ind_dup_no_singleton.RData")
save(type_dup,file="type_dup_no_singleton.RData")
save(pathway_data,file="pathway_data_no_singleton.RData")
save(ER,file="ER.RData")

## save all data
#save(d_nodup, featureNames_nodup, type_nodup, d_dup, uni_genes, pathway_ind_dup, fake_gene_ind_dup,
#	true_gene_ind_dup, type_dup, pathway_data, file="paper_app1_8pathways_input.RData")
