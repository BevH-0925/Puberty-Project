#Setup----
library(pheatmap) 
library(tidyverse)
library(ggplotify) 
library(heatmaply)
library(gplots)
library("viridis")
library(dendextend)

setwd("~/Documents/Code/RNASeqBH1-BH32/Random Forest")

#Load and prep data----
#Updated clinical data with new IDs & summary diagnosis
clindata<- read.csv("~/Documents/Code/RMCH Clinical Data/Data/Clindata_T.csv")
rownames(clindata) <- clindata$ParticipantNumber
IDs <- read.csv("~/Documents/Code/RMCH Clinical Data/Data/IDs.csv")
TIDs <- IDs[c(1:32),-1]

#Old ID column added to new clinical data
clindata <- cbind(TIDs,clindata)

#Load Prediction Accuracy column
MFglm_pred<- read.csv("~/Documents/Code/RNASeqBH1-BH32/Random Forest/MF_glm/MFglmIx_predictionaccuracy.csv")
#Put patients in same order as clinical data file
MFglm_pred <- MFglm_pred[match(clindata$TranscriptomeID, MFglm_pred$SampleID), ]
MFglm_pred <- MFglm_pred[,c(2,14)]
#Add prediction column to clinical data
clindata <- cbind(clindata,MFglm_pred)

#Load gene names from glm model
MFglm<- read.csv("~/Documents/Code/RNASeqBH1-BH32/Random Forest/MFglm_borcombined.csv")

Exp <- read.csv("~/Documents/Code/RNASeqBH1-BH32/EdgeR/Data/BH1-32T_lcpm.csv",row.names=1)
Exp <- t(Exp)
Exp <- as.data.frame(Exp)
Exp$rowID <- rownames(Exp)
#Put patients in same order as TID file
Exp <- Exp[match(TIDs$TranscriptomeID, Exp$rowID), ]
Exp <- cbind(TIDs,Exp)
rownames(Exp) <- Exp$NewIDforthesis
Exp <- Exp[,c(-2,-19588)]

MFglmExp <- Exp[,match(MFglm$Gene,colnames(Exp))]

clindata$ID<- paste0("(",clindata$ParticipantNumber, ")")
clindata$Xlabel <- paste(clindata$Heatmap_summary,clindata$ID)
tMFglmExp <- t(MFglmExp)
colnames(tMFglmExp) <- clindata$Xlabel

colnames(clindata)[colnames(clindata) == "Group"] <- "Pubertal Status"
colnames(clindata)[colnames(clindata) == "PredictionAccuracy"] <- "Prediction Accuracy"

cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

tMFglmExp_norm <- t(apply(tMFglmExp, 1, cal_z_score))
my_hclust_gene <- hclust(dist(tMFglmExp_norm), method = "complete")
my_gene_col <- cutree(tree = as.dendrogram(my_hclust_gene), k = 2)
my_gene_col <- data.frame(Cluster = ifelse(test = my_gene_col == 1, yes = "Cluster 1", no = "Cluster 2"))
my_sample_col <- clindata[,c(6,5,35,31,37)]
rownames(my_sample_col) <- my_sample_col$Xlabel
my_sample_col <- my_sample_col[,-5]
MFglmHM <- pheatmap(tMFglmExp_norm, annotation_row = my_gene_col, annotation_col = my_sample_col,  cutree_rows = 2,
                    cutree_cols = 5, fontsize_row = 4, fontsize_col = 5)

save_pheatmap_png <- function(x, filename, width=2400, height=2000, res = 300) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

save_pheatmap_png(MFglmHM, "MFglmHM.png")
