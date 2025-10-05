setwd("~/Documents/ALSPAC_Hypergraphs")

library(BioQC)
library(gplots)
library("dendextend") 

#Females closest to P2 (outliers removed)----
Top500_P2FCl <- read.csv("~/Documents/ALSPAC_DGE/Data/Expdata_FCl_or.csv")
rownames(Top500_P2FCl)<-paste(Top500_P2FCl[,1])
Top500_P2FCl <- Top500_P2FCl[,-1]

load("~/Documents/Code/EXP.Rdata")

#Transpose transcriptomic data so individuals are in rows, genes in columns
Expdata_t <- (t(data))
rm(data)

#Remove row 642 and 728 due to -9; Patient IDs #----A and ----A
Expdata_t <- Expdata_t[-c(642,728),]

#Subset of transcriptomic data to include only those participants present in Top500_P2FCl_or file
All_P2FCl <- Expdata_t[match(rownames(Top500_P2FCl),rownames(Expdata_t)),]

#Check both files in same order
rownames(All_P2FCl)[c(10,65,101,175,248)]
rownames(Top500_P2FCl)[c(10,65,101,175,248)]

#Positions of columns for top 500 probes in 'all probes' object
Positions <- match(colnames(Top500_P2FCl),colnames(All_P2FCl))

#All probes minus DEGs
MinusTop500_P2FCl <- All_P2FCl[,(-Positions)]

#check
which(colnames(MinusTop500_P2FCl)== "ILMN_1762337")
which(colnames(All_P2FCl)== "ILMN_1762337")
which(colnames(Top500_P2FCl)== "ILMN_1762337")

#Correlate DEGs against all other genes
cor_data_P2FCl<-cor(Top500_P2FCl,MinusTop500_P2FCl)

#Binarise correlation values using sd of correlation matrix
bin_P2FCl<-abs(cor_data_P2FCl)#remove sign
bin_P2FCl[which(bin_P2FCl>sd(cor_data_P2FCl))]<-1
bin_P2FCl[which(bin_P2FCl!=1)]<-0 #incidence matrix

#Matrix multiplication to generate hypernetwork adjacency matrix
hyp_P2FCl<-bin_P2FCl %*% t(bin_P2FCl) #adjacency matrix

#Generate heatmap
hm_P2FCl<-heatmap.2(hyp_P2FCl,trace="none", main="Females Closest to P2", labRow = FALSE, labCol=FALSE)

#ID genes in the central cluster
dend_P2FCl<-as.hclust(hm_P2FCl$rowDendrogram)#converts dendrogram part of hm object to class hclust
k_P2FCl<-2 #specify number of clusters based on dendrogram
ct_P2FCl<- cutree(dend_P2FCl, k_P2FCl)#cuts tree into groups, object tells you which cluster each probe is in
plot(dend_P2FCl,cex=0.5) 
rect.hclust(dend_P2FCl,k_P2FCl,which=1) #Visualize clusters on dendrogram, draws box around cluster x
length(ct_P2FCl[which(ct_P2FCl==2)]) #check cluster length is as expected for cluster x
clusterx_P2FCl<-names(ct_P2FCl[which(ct_P2FCl==1)]) #extract gene names for the selected cluster
write.csv(clusterx_P2FCl,file="Cluster1of2_P2FCl.csv") #csv file of gene names for the selected cluster

#Denextend - colours
P2FCl_dend_obj <- as.dendrogram(dend_P2FCl)
P2FCl_dend_obj%>% set("branches_k_color", k = 2) %>% plot(main = "Females Closest to P2 - into 2 clusters")
P2FCl_dend_obj %>% rect.dendrogram(k=2, border = 8, lty = 5, lwd = 2)

#Males closest to P2 (outliers removed)----
Top500_P2MCl <- read.csv("~/Documents/ALSPAC_DGE/Data/Expdata_MCl_or.csv")
rownames(Top500_P2MCl)<-paste(Top500_P2MCl[,1])
Top500_P2MCl <- Top500_P2MCl[,-1]

#Subset of transcriptomic data to include only those participants present in Top500_P2MCl_or file
All_P2MCl <- Expdata_t[match(rownames(Top500_P2MCl),rownames(Expdata_t)),]

#Check both files in same order
rownames(All_P2MCl)[c(10,65,101,175,204)]
rownames(Top500_P2MCl)[c(10,65,101,175,204)]

#Positions of columns for top 500 probes in 'all probes' object
Positions <- match(colnames(Top500_P2MCl),colnames(All_P2MCl))

#All probes minus DEGs
MinusTop500_P2MCl <- All_P2MCl[,(-Positions)]

#Correlate DEGs against all other genes
cor_data_P2MCl<-cor(Top500_P2MCl,MinusTop500_P2MCl)

#Binarise correlation values using sd of correlation matrix
bin_P2MCl<-abs(cor_data_P2MCl)#remove sign
bin_P2MCl[which(bin_P2MCl>sd(cor_data_P2MCl))]<-1
bin_P2MCl[which(bin_P2MCl!=1)]<-0 #incidence matrix

#Matrix multiplication to generate hypernetwork adjacency matrix
hyp_P2MCl<-bin_P2MCl %*% t(bin_P2MCl) #adjacency matrix

#Generate heatmap
hm_P2MCl<-heatmap.2(hyp_P2MCl,trace="none", main="Males Closest to P2", labRow = FALSE, labCol=FALSE)

#ID genes in the central cluster
dend_P2MCl<-as.hclust(hm_P2MCl$rowDendrogram)
k_P2MCl<-2 
ct_P2MCl<- cutree(dend_P2MCl, k_P2MCl)
plot(dend_P2MCl,cex=0.5) 
rect.hclust(dend_P2MCl,k_P2MCl,which=2) 
length(ct_P2MCl[which(ct_P2MCl==2)]) 
clusterx_P2MCl<-names(ct_P2MCl[which(ct_P2MCl==1)]) 
write.csv(clusterx_P2MCl,file="Cluster1of2_P2MCl.csv")

#Denextend - colours
P2MCl_dend_obj <- as.dendrogram(dend_P2MCl)
P2MCl_dend_obj%>% set("branches_k_color", k = 2) %>% plot(main = "Males Closest to P2 - into 2 clusters")
P2MCl_dend_obj %>% rect.dendrogram(k=2, border = 8, lty = 5, lwd = 2)

#Females furthest from P2 (outliers removed)----
Top500_P2FFu <- read.csv("~/Documents/ALSPAC_DGE/Data/Expdata_FFu_or.csv")
rownames(Top500_P2FFu)<-paste(Top500_P2FFu[,1])
Top500_P2FFu <- Top500_P2FFu[,-1]

#Subset of transcriptomic data to include only those participants present in Top500_P2FFu_or file
All_P2FFu <- Expdata_t[match(rownames(Top500_P2FFu),rownames(Expdata_t)),]

#Check both files in same order
rownames(All_P2FFu)[c(10,65,101,175,204)]
rownames(Top500_P2FFu)[c(10,65,101,175,204)]

#Positions of columns for top 500 probes in 'all probes' object
Positions <- match(colnames(Top500_P2FFu),colnames(All_P2FFu))

#All probes minus DEGs
MinusTop500_P2FFu <- All_P2FFu[,(-Positions)]

#Correlate DEGs against all other genes
cor_data_P2FFu<-cor(Top500_P2FFu,MinusTop500_P2FFu)

#Binarise correlation values using sd of correlation matrix
bin_P2FFu<-abs(cor_data_P2FFu)#remove sign
bin_P2FFu[which(bin_P2FFu>sd(cor_data_P2FFu))]<-1
bin_P2FFu[which(bin_P2FFu!=1)]<-0 #incidence matrix

#Matrix multiplication to generate hypernetwork adjacency matrix
hyp_P2FFu<-bin_P2FFu %*% t(bin_P2FFu) #adjacency matrix

#Generate heatmap
hm_P2FFu<-heatmap.2(hyp_P2FFu,trace="none", main="Females Furthest from P2", labRow = FALSE, labCol=FALSE)

#ID genes in the central cluster
dend_P2FFu<-as.hclust(hm_P2FFu$rowDendrogram)
k_P2FFu<-2 
ct_P2FFu<- cutree(dend_P2FFu, k_P2FFu)
plot(dend_P2FFu,cex=0.5) 
rect.hclust(dend_P2FFu,k_P2FFu,which=2)
length(ct_P2FFu[which(ct_P2FFu==2)]) 
clusterx_P2FFu<-names(ct_P2FFu[which(ct_P2FFu==2)]) 
write.csv(clusterx_P2FFu,file="Cluster2of2_P2FFu.csv")

#Denextend - colours
P2FFu_dend_obj <- as.dendrogram(dend_P2FFu)
P2FFu_dend_obj%>% set("branches_k_color", k = 2) %>% plot(main = "Females Furthest from P2 - into 2 clusters")
P2FFu_dend_obj %>% rect.dendrogram(k=2, border = 8, lty = 5, lwd = 2)

#Males furthest from P2 (outliers removed)----
Top500_P2MFu <- read.csv("~/Documents/ALSPAC_DGE/Data/Expdata_MFu_or.csv")
rownames(Top500_P2MFu)<-paste(Top500_P2MFu[,1])
Top500_P2MFu <- Top500_P2MFu[,-1]

#Subset of transcriptomic data to include only those participants present in Top500_P2MFu_or file
All_P2MFu <- Expdata_t[match(rownames(Top500_P2MFu),rownames(Expdata_t)),]

#Check both files in same order
rownames(All_P2MFu)[c(10,65,101,175,204)]
rownames(Top500_P2MFu)[c(10,65,101,175,204)]

#Positions of columns for top 500 probes in 'all probes' object
Positions <- match(colnames(Top500_P2MFu),colnames(All_P2MFu))

#All probes minus DEGs
MinusTop500_P2MFu <- All_P2MFu[,(-Positions)]

#Correlate DEGs against all other genes
cor_data_P2MFu<-cor(Top500_P2MFu,MinusTop500_P2MFu)

#Binarise correlation values using sd of correlation matrix
bin_P2MFu<-abs(cor_data_P2MFu)#remove sign
bin_P2MFu[which(bin_P2MFu>sd(cor_data_P2MFu))]<-1
bin_P2MFu[which(bin_P2MFu!=1)]<-0 #incidence matrix

#Matrix multiplication to generate hypernetwork adjacency matrix
hyp_P2MFu<-bin_P2MFu %*% t(bin_P2MFu) #adjacency matrix

#Generate heatmap
hm_P2MFu<-heatmap.2(hyp_P2MFu,trace="none", main="Males Furthest from P2", labRow = FALSE, labCol=FALSE)

#ID genes in the central cluster
dend_P2MFu<-as.hclust(hm_P2MFu$rowDendrogram)
k_P2MFu<-2 #amend as required
ct_P2MFu<- cutree(dend_P2MFu, k_P2MFu)
plot(dend_P2MFu,cex=0.5) 
rect.hclust(dend_P2MFu,k_P2MFu,which=2) #amend as required
length(ct_P2MFu[which(ct_P2MFu==2)]) #amend as required
clusterx_P2MFu<-names(ct_P2MFu[which(ct_P2MFu==2)]) #amend as required
write.csv(clusterx_P2MFu,file="Cluster2of2_P2MFu.csv")

#Denextend - colours
P2MFu_dend_obj <- as.dendrogram(dend_P2MFu)
P2MFu_dend_obj%>% set("branches_k_color", k = 2) %>% plot(main = "Males Furthest from P2 - into 2 clusters") #amend as required
P2MFu_dend_obj %>% rect.dendrogram(k=2, border = 8, lty = 5, lwd = 2)