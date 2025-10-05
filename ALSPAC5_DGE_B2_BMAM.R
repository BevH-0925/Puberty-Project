setwd("~/Documents/Code")

library(mixOmics)
library(rgl)
library("limma")
library("Biobase")
library(tidyverse)
library(ggpubr)

load("~/Documents/Files from PEC drive/ALSPAC_Original Files/EXP.Rdata")

#Transpose transcriptomic data so individuals are in rows, genes in columns
Expdata_t <- (t(data))
rm(data)

#Remove row 642 and 728 due to -9; Patient IDs #----A and ----A
Expdata_t <- Expdata_t[-c(642,728),]

#Load clinical data with distance to B2 quartiles 
ClinData_B2_Quartiles <- read.csv("~/Documents/Code/B2quartilesT.csv")

#Remove column X.1 & add row names
ClinData_B2_Quartiles <- ClinData_B2_Quartiles[,-1]
rownames(ClinData_B2_Quartiles)<-paste(ClinData_B2_Quartiles[,1])
ClinData_B2_Quartiles <- ClinData_B2_Quartiles[,-1]

#Subset of transcriptomic data to include only those with a value for TrB2
Expdata_DtB2 <- Expdata_t[match(rownames(ClinData_B2_Quartiles),rownames(Expdata_t)),]
rm(Expdata_t)

#Check both files in same order
rownames(Expdata_DtB2)[c(10,65,275,361)]
rownames(ClinData_B2_Quartiles)[c(10,65,275,361)]

#Groups
#Time to B2 below median vs above median
MedianB2<- ClinData_B2_Quartiles[,c(1,2292)]
DtB2Median <- median(MedianB2$TrB2)
MedianB2$Group <- ifelse(MedianB2$TrB2<DtB2Median, "Below median", NA)
MedianB2$Group <- ifelse(MedianB2$TrB2>=DtB2Median, "Equal to or above median", MedianB2$Group)
MedianB2 <- MedianB2[,-1]

#Convert Group column to a factor
Y <- as.factor(MedianB2$Group)
summary(Y)

#PCA all probes, time to B2 below B2 vs above median group
Y.pca <- pca(Expdata_DtB2, ncomp=364) 
plotIndiv(Y.pca, group=Y, title="PCA, transcriptome, all probes", legend=TRUE, 
          legend.title="Time to B2", size.xlabel=rel(1.5), size.ylabel = rel(1.5), size.axis=rel(1.5), 
          size.legend.title=rel(1.5))

#Calculate variance
eigs <- Y.pca$sdev^2
PCA1_var <- eigs[1] / sum(eigs)
PCA1_var<-round(PCA1_var,digits = 2)*100
PCA2_var <- eigs[2] / sum(eigs)
PCA2_var<-round(PCA2_var,digits = 2)*100

pcadata_plot <- data.frame(Y.pca$variates$X[,1:2])

ggplot(pcadata_plot,aes(PC1,PC2,colour=Y, shape=Y))+
  geom_point(size=5, alpha=0.8)+
  coord_fixed(ratio=1)+
  theme_classic(base_size=30)+
  labs(title= "a) PCA, all probes", colour="Time to B2", shape="Time to B2")+
  xlab(label= paste("PC1 (expl. variance= ", PCA1_var, "%)",sep=""))+
  ylab(label= paste("PC2 (expl. variance= ", PCA2_var, "%)",sep=""))+
  theme(legend.position = "bottom")+
  scale_color_manual(values = c("#1F77B4", "#FF7F0E"), labels=c("BM", "AM"))+
  scale_shape_manual(values=c(1, 2), labels=c("BM", "AM"))
ggsave("SQaPCAallBMAM.png", width = 10, height = 10, dpi=320)

#OUTLIERS
#New clinical data files with outlier participants removed
#Removing outliers identified on original PCA
#Identify positions of outliers in clin data file
which(rownames(ClinData_B2_Quartiles) == "----A")
which(rownames(ClinData_B2_Quartiles) == "----A")
#etc.

#Create file with outliers removed
ClinData_B2_Quartiles_or <- ClinData_B2_Quartiles[-c(220, 267, 310, 309, 303, 192, 232, 178, 136, 208, 258, 343, 261, 285, 274, 129, 337, 177, 10, 211, 2, 157,266,315,320),]

#Recalculation of median now that outliers removed
MedianB2_or<- ClinData_B2_Quartiles_or[,c(1,2292)]
DtB2Median_or <- median(MedianB2_or$TrB2)
MedianB2_or$Group <- ifelse(MedianB2_or$TrB2<DtB2Median_or, "Below median", NA)
MedianB2_or$Group <- ifelse(MedianB2_or$TrB2>=DtB2Median_or, "Equal to or above median", MedianB2_or$Group)
MedianB2_or <- MedianB2_or[,-1]

#Convert Group column to a factor
Y_or <- as.factor(MedianB2_or$Group)
summary(Y_or) #below median 161, equal to or above median 178

#Subset of transcriptomic data to remove outliers
Expdata_DtB2_or <- Expdata_DtB2[match(rownames(ClinData_B2_Quartiles_or),rownames(Expdata_DtB2)),]

#PLSDA, all probes, time to B2 below median vs above median group, outliers removed
plsda.data<- plsda(Expdata_DtB2_or, Y_or, ncomp=339)

#Calculate variance
ev <- explained_variance(plsda.data$X, plsda.data$variates$X, ncomp=339)
Comp1v<-round(ev[1],digits = 2)*100
Comp2v<-round(ev[2],digits = 2)*100

plsdadata_plot <- data.frame(plsda.data$variates$X[,1:2])

ggplot(plsdadata_plot,aes(comp1,comp2,colour=Y_or, shape=Y_or))+
  geom_point(size=5, alpha=0.8)+
  coord_fixed(ratio=1)+
  theme_classic(base_size=30)+
  labs(title= "d) PLSDA, all probes", colour="Time to B2", shape="Time to B2")+
  xlab(label= paste("Component 1 (expl. variance= ", Comp1v, "%)",sep=""))+
  ylab(label= paste("Component 2 (expl. variance= ", Comp2v, "%)",sep=""))+
  theme(legend.position= "bottom")+
  scale_color_manual(values = c("#1F77B4", "#FF7F0E"), labels=c("BM", "AM"))+
  scale_shape_manual(values=c(1, 2), labels=c("BM", "AM"))+
  stat_stars(aes(colour=Y_or), alpha=0.5)+
  stat_ellipse((aes(colour=Y_or)), level=0.95, linewidth = 1.5, type="t", alpha=0.5)
ggsave("SQdPLSDAallBMAM.png", width = 10, height = 10, dpi=320)

#PCA all probes, time to B2 below median vs above median group, outliers removed
Y_or.pca <- pca(Expdata_DtB2_or, ncomp=339)

#Calculate variance
eigs <- Y_or.pca $sdev^2
PCA1_var <- eigs[1] / sum(eigs)
PCA1_var<-round(PCA1_var,digits = 2)*100
PCA2_var <- eigs[2] / sum(eigs)
PCA2_var<-round(PCA2_var,digits = 2)*100

pcadata_plot <- data.frame(Y_or.pca$variates$X[,1:2])

ggplot(pcadata_plot,aes(PC1,PC2,colour=Y_or, shape=Y_or))+
  geom_point(size=5, alpha=0.8)+
  coord_fixed(ratio=1)+
  theme_classic(base_size=30)+
  labs(title= "b) PCA, all probes, outliers removed", colour="Time to B2", shape="Time to B2")+
  xlab(label= paste("PC1 (expl. variance= ", PCA1_var, "%)",sep=""))+
  ylab(label= paste("PC2 (expl. variance= ", PCA2_var, "%)",sep=""))+
  theme(legend.position = "bottom")+
  scale_color_manual(values = c("#1F77B4", "#FF7F0E"), labels=c("BM", "AM"))+
  scale_shape_manual(values=c(1, 2), labels=c("BM", "AM"))+
  stat_stars(aes(colour=Y_or), alpha=0.5)+
  stat_ellipse((aes(colour=Y_or)), level=0.95, linewidth = 1.5, type="t", alpha=0.5)
ggsave("SQbPCAallorBMAM.png", width = 10, height = 10, dpi=320)

#Expression set - time to B2 below median vs above median group (outliers removed)
Eset_expdata_DtB2Median_or <- ExpressionSet(t(Expdata_DtB2_or))
phenoData(Eset_expdata_DtB2Median_or) <- AnnotatedDataFrame(MedianB2_or)
design_MedianB2_or <- model.matrix(~Eset_expdata_DtB2Median_or@phenoData$Group)
fit_MedianB2_or <- lmFit(Eset_expdata_DtB2Median_or,design_MedianB2_or)
fit_MedianB2_or <- eBayes(fit_MedianB2_or)
list_MedianB2_or <- topTable(fit_MedianB2_or, number = 500, adjust.method = "BH")
write.csv(list_MedianB2_or,file="Top500_MedianB2_or.csv")

#Analysis using top 500 differentially expressed probes
#Create vector of probe names for each group
Probes_MedianB2_or <- rownames(list_MedianB2_or)

#Create subset of transcription data with top 500 differentially expressed probes and outlier participants removed
Expdata_MedianB2_or <- Expdata_DtB2_or[,match(Probes_MedianB2_or,colnames(Expdata_DtB2_or))]

write.csv(Expdata_MedianB2_or,file="Expdata_MedianB2_or.csv")

#PLSDA top 500 DGE, time to B2 below median vs above median group (outliers removed)
plsda.dataY_or<- plsda(Expdata_MedianB2_or, Y_or, ncomp = 339)

#Calculate variance
ev <- explained_variance(plsda.dataY_or$X, plsda.dataY_or$variates$X, ncomp=339)
Comp1v<-round(ev[1],digits = 2)*100
Comp2v<-round(ev[2],digits = 2)*100

plsdadata_plot <- data.frame(plsda.dataY_or$variates$X[,1:2])

ggplot(plsdadata_plot,aes(comp1,comp2,colour=Y_or, shape=Y_or))+
  geom_point(size=5, alpha=0.8)+
  coord_fixed(ratio=1)+
  theme_classic(base_size=30)+
  labs(title= "e) PLSDA, top 500", colour="Time to B2", shape="Time to B2")+
  xlab(label= paste("Component 1 (expl. variance= ", Comp1v, "%)",sep=""))+
  ylab(label= paste("Component 2 (expl. variance= ", Comp2v, "%)",sep=""))+
  theme(legend.position= "bottom")+
  scale_color_manual(values = c("#1F77B4", "#FF7F0E"), labels=c("BM", "AM"))+
  scale_shape_manual(values=c(1, 2), labels=c("BM", "AM"))+
  stat_stars(aes(colour=Y_or), alpha=0.5)+
  stat_ellipse((aes(colour=Y_or)), level=0.95, linewidth = 1.5, type="t", alpha=0.5)
ggsave("SQePLSDA500BMAM.png", width = 10, height = 10, dpi=320)

#PCA top 500 DGE, time to B2 below median vs above median group (outliers removed)
Y_or.pca <- pca(Expdata_MedianB2_or, ncomp=339) 

#Calculate variance
eigs <- Y_or.pca $sdev^2
PCA1_var <- eigs[1] / sum(eigs)
PCA1_var<-round(PCA1_var,digits = 2)*100
PCA2_var <- eigs[2] / sum(eigs)
PCA2_var<-round(PCA2_var,digits = 2)*100

pcadata_plot <- data.frame(Y_or.pca$variates$X[,1:2])

ggplot(pcadata_plot,aes(PC1,PC2,colour=Y_or, shape=Y_or))+
  geom_point(size=5, alpha=0.8)+
  coord_fixed(ratio=1)+
  theme_classic(base_size=30)+
  labs(title= "c) PCA, top 500", colour="Time to B2", shape="Time to B2")+
  xlab(label= paste("PC1 (expl. variance= ", PCA1_var, "%)",sep=""))+
  ylab(label= paste("PC2 (expl. variance= ", PCA2_var, "%)",sep=""))+
  theme(legend.position = "bottom")+
  scale_color_manual(values = c("#1F77B4", "#FF7F0E"), labels=c("BM", "AM"))+
  scale_shape_manual(values=c(1, 2), labels=c("BM", "AM"))+
  stat_stars(aes(colour=Y_or), alpha=0.5)+
  stat_ellipse((aes(colour=Y_or)), level=0.95, linewidth = 1.5, type="t", alpha=0.5)
ggsave("SQcPCA500BMAM.png", width = 10, height = 10, dpi=320)

write.csv(MedianB2_or,file="MedianB2_or.csv")
