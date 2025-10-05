library(BioQC)
library(gplots)
library("dendextend")
library(ggplot2)

setwd("~/Documents/RNASeqBH1-BH32/Hypergraphs")

#All 32 participants classic EdgeR list----
Top500_MFcl <- read.csv("~/Documents/RNASeqBH1-BH32/EdgeR/Data/BH1-32T_top500_lcpm_classic.csv")
rownames(Top500_MFcl)<-paste(Top500_MFcl[,1])
Top500_MFcl <- Top500_MFcl[,-1]
Top500_MFcl <- (t(Top500_MFcl))

lcpm_all <- read.csv("~/Documents/RNASeqBH1-BH32/EdgeR/Data/BH1-32T_lcpm.csv")
rownames(lcpm_all)<-paste(lcpm_all[,1])
lcpm_all<- lcpm_all[,-1]
lcpm_all <- (t(lcpm_all))

#Positions of columns for top 500 probes in 'all probes' object
Positions <- match(colnames(Top500_MFcl),colnames(lcpm_all))

#All probes minus DEGs
MinusTop500_MFcl <- lcpm_all[,(-Positions)]

#check
which(colnames(MinusTop500_MFcl)== "DPM1")
which(colnames(lcpm_all)== "DPM1")
which(colnames(Top500_MFcl)== "DPM1")

#Correlate DEGs against all other genes
cor_data_MFcl<-cor(Top500_MFcl,MinusTop500_MFcl)

#Binarise correlation values using sd of correlation matrix
bin_MFcl<-abs(cor_data_MFcl)#remove sign
bin_MFcl[which(bin_MFcl>sd(cor_data_MFcl))]<-1
bin_MFcl[which(bin_MFcl!=1)]<-0 #incidence matrix

#Matrix multiplication to generate hypernetwork adjacency matrix
hyp_MFcl<-bin_MFcl %*% t(bin_MFcl) #adjacency matrix

#Generate heatmap
hm_MFcl<-heatmap.2(hyp_MFcl,trace="none", labRow = FALSE, labCol=FALSE)

#ID genes in the central cluster
dend_MFcl<-as.hclust(hm_MFcl$rowDendrogram)#converts dendrogram part of hm object to class hclust
k_MFcl<-2 #specify number of clusters based on dendrogram
ct_MFcl<- cutree(dend_MFcl, k_MFcl)#cuts tree into groups, object tells you which cluster each probe is in
plot(dend_MFcl,cex=0.5) #plotted dendrogram can be misleading, ID of clusters can differ between ct and the dend plot
rect.hclust(dend_MFcl,k_MFcl,which=2) #amend as required
length(ct_MFcl[which(ct_MFcl==1)]) #amend as required
clusterx_MFcl<-names(ct_MFcl[which(ct_MFcl==2)]) #amend as required
write.csv(clusterx_MFcl,file="Cluster1of2_MFcl.csv") 

#Calculate hypergraph row sums
hgrs <- rowSums(hyp_MFcl)
write.csv(hgrs,file="MFcl_hypergraph_rowsum.csv")

#Denextend - colours
MFcl_dend_obj <- as.dendrogram(dend_MFcl)
MFcl_dend_obj%>% set("branches_k_color", k = 2) %>% plot(main = "BH1-32 MF classic - into 2 clusters")
MFcl_dend_obj %>% rect.dendrogram(k=2, border = 8, lty = 5, lwd = 2)

#Female participants classic EdgeR list----
Top500_Fcl <- read.csv("~/Documents/RNASeqBH1-BH32/EdgeR/Data/BH1-32T_top500_lcpm_classic_F.csv")
rownames(Top500_Fcl)<-paste(Top500_Fcl[,1])
Top500_Fcl <- Top500_Fcl[,-1]
Top500_Fcl<- (t(Top500_Fcl))

lcpm_F <- read.csv("~/Documents/RNASeqBH1-BH32/EdgeR/Data/BH1-32T_lcpm_F.csv")
rownames(lcpm_F)<-paste(lcpm_F[,1])
lcpm_F<- lcpm_F[,-1]
lcpm_F <- (t(lcpm_F))

#Subset of transcriptomic data to include only those participants present in Top500 object
lcpm_Fsub <- lcpm_F[match(rownames(Top500_Fcl),rownames(lcpm_F)),]

#Check both files in same order
rownames(lcpm_Fsub)[c(2,5,11,20)]
rownames(Top500_Fcl)[c(2,5,11,20)]

#Positions of columns for top 500 probes in 'all probes' object
Positions <- match(colnames(Top500_Fcl),colnames(lcpm_Fsub))

#All probes minus DEGs
MinusTop500_Fcl <- lcpm_Fsub[,(-Positions)]

#check
which(colnames(MinusTop500_Fcl)== "PPP1R3G")
which(colnames(lcpm_Fsub)== "PPP1R3G")
which(colnames(Top500_Fcl)== "PPP1R3G")

#Correlate DEGs against all other genes
cor_data_Fcl<-cor(Top500_Fcl,MinusTop500_Fcl)

#Binarise correlation values using sd of correlation matrix
bin_Fcl<-abs(cor_data_Fcl)#remove sign
bin_Fcl[which(bin_Fcl>sd(cor_data_Fcl))]<-1
bin_Fcl[which(bin_Fcl!=1)]<-0 #incidence matrix

#Matrix multiplication to generate hypernetwork adjacency matrix
hyp_Fcl<-bin_Fcl %*% t(bin_Fcl) #adjacency matrix

ent_Fcl<-entropy(hyp_Fcl)

#Generate heatmap
hm_Fcl<-heatmap.2(hyp_Fcl,trace="none", labRow = FALSE, labCol=FALSE)

#ID genes in the central cluster
dend_Fcl<-as.hclust(hm_Fcl$rowDendrogram)#converts dendrogram part of hm object to class hclust
k_Fcl<-4 #specify number of clusters based on dendrogram
ct_Fcl<- cutree(dend_Fcl, k_Fcl)#cuts tree into groups, object tells you which cluster each probe is in
plot(dend_Fcl,cex=0.5) #plotted dendrogram can be misleading, ID of clusters can differ between ct and the dend plot
rect.hclust(dend_Fcl,k_Fcl,which=1) #Visualize clusters on dendrogram, draws box around cluster x
length(ct_Fcl[which(ct_Fcl==1)]) #check cluster length is as expected for cluster x
clusterx_Fcl<-names(ct_Fcl[which(ct_Fcl==4)]) #extract gene names for the selected cluster
write.csv(clusterx_Fcl,file="Cluster4of4_Fcl.csv") #csv file of gene names for the selected cluster

#Calculate hypergraph row sums
hgrs <- rowSums(hyp_Fcl)
write.csv(hgrs,file="Fcl_hypergraph_rowsum.csv")

#Denextend - colours
Fcl_dend_obj <- as.dendrogram(dend_Fcl)
Fcl_dend_obj%>% set("branches_k_color", k = 2) %>% plot(main = "BH1-32 F classic - into 2 clusters")
Fcl_dend_obj %>% rect.dendrogram(k=2, border = 8, lty = 5, lwd = 2)

#All 32 participants glm adjusted EdgeR list----
Top500_MFglm<- read.csv("~/Documents/RNASeqBH1-BH32/EdgeR/Data/BH1-32T_top500_lcpm_glmsea.csv")
rownames(Top500_MFglm)<-paste(Top500_MFglm[,1])
Top500_MFglm <- Top500_MFglm[,-1]
Top500_MFglm <- (t(Top500_MFglm))

#Positions of columns for top 500 probes in 'all probes' object
Positions <- match(colnames(Top500_MFglm),colnames(lcpm_all))

#All probes minus DEGs
MinusTop500_MFglm <- lcpm_all[,(-Positions)]

#check
which(colnames(MinusTop500_MFglm)== "DPM1")
which(colnames(lcpm_all)== "DPM1")
which(colnames(Top500_MFglm)== "DPM1")

#Correlate DEGs against all other genes
cor_data_MFglm<-cor(Top500_MFglm,MinusTop500_MFglm)

#Binarise correlation values using sd of correlation matrix
bin_MFglm<-abs(cor_data_MFglm)#remove sign
bin_MFglm[which(bin_MFglm>sd(cor_data_MFglm))]<-1
bin_MFglm[which(bin_MFglm!=1)]<-0 #incidence matrix

#Matrix multiplication to generate hypernetwork adjacency matrix
hyp_MFglm<-bin_MFglm %*% t(bin_MFglm) #adjacency matrix

#Generate heatmap
hm_MFglm<-heatmap.2(hyp_MFglm,trace="none", labRow = FALSE, labCol=FALSE)

#ID genes in the central cluster
dend_MFglm<-as.hclust(hm_MFglm$rowDendrogram)#converts dendrogram part of hm object to class hclust
k_MFglm<-2 #specify number of clusters based on dendrogram
ct_MFglm<- cutree(dend_MFglm, k_MFglm)#cuts tree into groups, object tells you which cluster each probe is in
plot(dend_MFglm,cex=0.5) #plotted dendrogram can be misleading, ID of clusters can differ between ct and the dend plot
rect.hclust(dend_MFglm,k_MFglm,which=2) #Visualize clusters on dendrogram, draws box around cluster x
length(ct_MFglm[which(ct_MFglm==1)]) #check cluster length is as expected for cluster x
clusterx_MFglm<-names(ct_MFglm[which(ct_MFglm==1)]) #extract gene names for the selected cluster
write.csv(clusterx_MFglm,file="Cluster1of2_MFglm.csv") #csv file of gene names for the selected cluster

#Calculate hypergraph row sums
hgrs <- rowSums(hyp_MFglm)
write.csv(hgrs,file="MFglm_hypergraph_rowsum.csv")

#Denextend - colours
MFglm_dend_obj <- as.dendrogram(dend_MFglm)
MFglm_dend_obj%>% set("branches_k_color", k = 2) %>% plot(main = "BH1-32 MF glm - into 2 clusters")
MFglm_dend_obj %>% rect.dendrogram(k=2, border = 8, lty = 5, lwd = 2)

#Female participants glm EdgeR list----
Top500_Fglm <- read.csv("~/Documents/RNASeqBH1-BH32/EdgeR/Data/BH1-32T_top500_lcpm_F_glmea.csv")
rownames(Top500_Fglm)<-paste(Top500_Fglm[,1])
Top500_Fglm <- Top500_Fglm[,-1]
Top500_Fglm<- (t(Top500_Fglm))

#Subset of transcriptomic data to include only those participants present in Top500 object
lcpm_Fsub <- lcpm_F[match(rownames(Top500_Fglm),rownames(lcpm_F)),]

#Check both files in same order
rownames(lcpm_Fsub)[c(2,5,11,20)]
rownames(Top500_Fglm)[c(2,5,11,20)]

#Positions of columns for top 500 probes in 'all probes' object
Positions <- match(colnames(Top500_Fglm),colnames(lcpm_Fsub))

#All probes minus DEGs
MinusTop500_Fglm <- lcpm_Fsub[,(-Positions)]

#check
which(colnames(MinusTop500_Fglm)== "PPP1R3G")
which(colnames(lcpm_Fsub)== "PPP1R3G")
which(colnames(Top500_Fglm)== "PPP1R3G")

#Correlate DEGs against all other genes
cor_data_Fglm<-cor(Top500_Fglm,MinusTop500_Fglm)

#Binarise correlation values using sd of correlation matrix
bin_Fglm<-abs(cor_data_Fglm)#remove sign
bin_Fglm[which(bin_Fglm>sd(cor_data_Fglm))]<-1
bin_Fglm[which(bin_Fglm!=1)]<-0 #incidence matrix

#Matrix multiplication to generate hypernetwork adjacency matrix
hyp_Fglm<-bin_Fglm %*% t(bin_Fglm) #adjacency matrix

#Generate heatmap
hm_Fglm<-heatmap.2(hyp_Fglm,trace="none", labRow = FALSE, labCol=FALSE)

#ID genes in the central cluster
dend_Fglm<-as.hclust(hm_Fglm$rowDendrogram)#converts dendrogram part of hm object to class hclust
k_Fglm<-2 #specify number of clusters based on dendrogram
ct_Fglm<- cutree(dend_Fglm, k_Fglm)#cuts tree into groups, object tells you which cluster each probe is in
plot(dend_Fglm,cex=0.5)
rect.hclust(dend_Fglm,k_Fglm,which=1) #Visualize clusters on dendrogram, draws box around cluster x
length(ct_Fglm[which(ct_Fglm==1)]) #check cluster length is as expected for cluster x
clusterx_Fglm<-names(ct_Fglm[which(ct_Fglm==1)]) #extract gene names for the selected cluster
write.csv(clusterx_Fglm,file="Cluster1of2_Fglm.csv") #csv file of gene names for the selected cluster

#Calculate hypergraph row sums
hgrs <- rowSums(hyp_Fglm)
write.csv(hgrs,file="Fglm_hypergraph_rowsum.csv")

#Denextend - colours
Fglm_dend_obj <- as.dendrogram(dend_Fglm)
Fglm_dend_obj%>% set("branches_k_color", k = 2) %>% plot(main = "BH1-32 F glm - into 2 clusters")
Fglm_dend_obj %>% rect.dendrogram(k=2, border = 8, lty = 5, lwd = 2)
