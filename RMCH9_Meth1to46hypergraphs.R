library(BioQC)
library(gplots)
library("dendextend")
library(ggplot2)

#All 46 methylome participants----
#Top 500 subset of normalised betas
Top500_MF <- MyNormt[,match(rownames(DMPtop500),colnames(MyNormt))]

#Positions of columns for top 500 probes in 'all probes' object
Positions <- match(colnames(Top500_MF),colnames(MyNormt))

#All probes minus DMPs
MinusTop500_MF <- MyNormt[,(-Positions)]

#check
which(colnames(MinusTop500_MF)== "cg17315310_BC21")
which(colnames(MyNormt)== "cg17315310_BC21")
which(colnames(Top500_MF)== "cg17315310_BC21")

#Correlate DEGs against all other genes
cor_data_MF<-cor(Top500_MF,MinusTop500_MF)

#Binarise correlation values using sd of correlation matrix
bin_MF<-abs(cor_data_MF)#remove sign
bin_MF[which(bin_MF>sd(cor_data_MF))]<-1
bin_MF[which(bin_MF!=1)]<-0 #incidence matrix

#Matrix multiplication to generate hypernetwork adjacency matrix
hyp_MF<-bin_MF %*% t(bin_MF) #adjacency matrix

#Generate heatmap
hm_MF<-heatmap.2(hyp_MF,trace="none", labRow = FALSE, labCol=FALSE)

#ID genes in the central cluster
dend_MF<-as.hclust(hm_MF$rowDendrogram)#converts dendrogram part of hm object to class hclust
k_MF<-3 #specify number of clusters based on dendrogram
ct_MF<- cutree(dend_MF, k_MF)#cuts tree into groups, object tells you which cluster each probe is in
plot(dend_MF,cex=0.5) #plotted dendrogram can be misleading, ID of clusters can differ between ct and the dend plot
rect.hclust(dend_MF,k_MF,which=3) #Visualize clusters on dendrogram, draws box around cluster x
length(ct_MF[which(ct_MF==3)]) #check cluster length is as expected for cluster x
clusterx_MF<-names(ct_MF[which(ct_MF==3)]) #extract gene names for the selected cluster
write.csv(clusterx_MF,file="Cluster3of3_MF.csv") #csv file of gene names for the selected cluster

#Calculate hypergraph row sums
hgrs <- rowSums(hyp_MF)
write.csv(hgrs,file="MF_hypergraph_rowsum.csv")
Entvsrowsum <- cbind(hgrs,Ent_hyp_MF,Entropy_rank)

ggplot(Entvsrowsum, aes(hgrs,Ent_hyp_MF))+
  geom_point()+
  labs(title="Batch 1 Methylome",x="Hypergraph Row Sum", y="Gene Entropy")+
  theme_classic()
ggsave("EntropyvsHGRS_MF.png")

ggplot(Entvsrowsum, aes(hgrs,Entropy_rank))+
  geom_point()+
  labs(title="Batch 1 Methylome",x="Hypergraph Row Sum", y="Gene Entropy Rank")+
  theme_classic()
ggsave("EntropyRankvsHGRS_MF.png")

#Denextend - colours
MF_dend_obj <- as.dendrogram(dend_MF)
MF_dend_obj%>% set("branches_k_color", k = 3) %>% plot(main = "Batch 1 Methylome - into 3 clusters")
MF_dend_obj %>% rect.dendrogram(k=3, border = 8, lty = 5, lwd = 2)

#32 F methylome participants----
#Top 500 subset of normalised betas
Top500_F <- MyNormFt[,match(DMPtop500F$X,colnames(MyNormFt))]

#Positions of columns for top 500 probes in 'all probes' object
Positions <- match(colnames(Top500_F),colnames(MyNormFt))

#All probes minus DMPs
MinusTop500_F <- MyNormFt[,(-Positions)]

#check
which(colnames(MinusTop500_F)== "cg15756399_TC21")
which(colnames(MyNormFt)== "cg15756399_TC21")
which(colnames(Top500_F)== "cg15756399_TC21")

#Correlate DEGs against all other genes
cor_data_F<-cor(Top500_F,MinusTop500_F)

#Binarise correlation values using sd of correlation matrix
bin_F<-abs(cor_data_F)#remove sign
bin_F[which(bin_F>sd(cor_data_F))]<-1
bin_F[which(bin_F!=1)]<-0 #incidence matrix

#Matrix multiplication to generate hypernetwork adjacency matrix
hyp_F<-bin_F %*% t(bin_F) #adjacency matrix

#Generate heatmap
hm_F<-heatmap.2(hyp_F,trace="none", labRow = FALSE, labCol=FALSE)

#ID genes in the central cluster
dend_F<-as.hclust(hm_F$rowDendrogram)#converts dendrogram part of hm object to class hclust
k_F<-4 #specify number of clusters based on dendrogram
ct_F<- cutree(dend_F, k_F)#cuts tree into groups, object tells you which cluster each probe is in
plot(dend_F,cex=0.5) #plotted dendrogram can be misleading, ID of clusters can differ between ct and the dend plot
rect.hclust(dend_F,k_F,which=1) #Visualize clusters on dendrogram, draws box around cluster x
length(ct_F[which(ct_F==1)]) #check cluster length is as expected for cluster x
clusterx_F<-names(ct_F[which(ct_F==4)]) #extract gene names for the selected cluster
write.csv(clusterx_F,file="Cluster4of4_F.csv") #csv file of gene names for the selected cluster

#Calculate hypergraph row sums
hgrs <- rowSums(hyp_F)
write.csv(hgrs,file="F_hypergraph_rowsum.csv")

#Denextend - colours
F_dend_obj <- as.dendrogram(dend_F)
F_dend_obj%>% set("branches_k_color", k = 4) %>% plot(main = "Batch 1 Methylome F- into 4 clusters")
F_dend_obj %>% rect.dendrogram(k=4, border = 8, lty = 5, lwd = 2)

