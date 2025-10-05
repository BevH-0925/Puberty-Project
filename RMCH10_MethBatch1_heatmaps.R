#Setup----
setwd("~/Documents/Code/RMCH_MethArrayBatch1/RMCHMethChAMP")

library(pheatmap) ## for heatmap generation
library(tidyverse) ## for data wrangling
library(ggplotify) ## to convert pheatmap to ggplot2
library(heatmaply) ## for constructing interactive heatmap
library(gplots)
library("viridis")
library(dendextend)

#Load and prep data----
#Updated clinical data with new IDs & summary diagnosis
clindata<- read.csv("~/Documents/Code/RMCH Clinical Data/Data/Clindata_Meth.csv")
rownames(clindata) <- clindata$ParticipantNumber
IDs <- read.csv("~/Documents/Code/RMCH Clinical Data/Data/IDs.csv")
MIDs <- IDs[,-2]

MIDs <-MIDs[match(clindata$ParticipantNumber, MIDs$NewIDforthesis ), ]
clindata <- cbind(MIDs,clindata)

MFIx_perf <- read.csv("~/Documents/Code/RMCH_MethArrayBatch1/RMCHMethChAMP/MFIx_perf.csv")
DMP_Pred <- read.csv("~/Documents/Code/RMCH_MethArrayBatch1/RMCHMethChAMP/DMP_Pred.csv")
load("~/Documents/Code/RMCH_MethArrayBatch1/RMCHMethChAMP/MyNormt.robj")

MFIx_perf <-MFIx_perf[,c(2,15)]
MFIx_perf <-MFIx_perf[match(clindata$ParticipantID, MFIx_perf$SampleID), ]
clindata <- cbind(MFIx_perf,clindata)
clindata$ID<- paste0("(",clindata$ParticipantNumber, ")")
clindata$Xlabel <- paste(clindata$Heatmap_summary,clindata$ID)

colnames(clindata)[colnames(clindata) == "Group"] <- "Pubertal Status"
colnames(clindata)[colnames(clindata) == "PredictionAccuracy"] <- "Prediction Accuracy"
rownames(clindata) <- clindata$SampleID

DMPPredbetas <- MyNormt[,match(DMP_Pred$Probe,colnames(MyNormt))]
DMPPredbetas$Xlabel <- clindata$Xlabel[match(rownames(DMPPredbetas), rownames(clindata))]
rownames(DMPPredbetas) <- DMPPredbetas$Xlabel
rownames(clindata) <- clindata$Xlabel
DMPPredbetas <- DMPPredbetas[,-171]

cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

DMPPredbetas_norm <- apply(DMPPredbetas, 1, cal_z_score)
my_hclust_gene <- hclust(dist(DMPPredbetas_norm), method = "complete")
my_gene_col <- cutree(tree = as.dendrogram(my_hclust_gene), k = 5)
my_gene_col<- as.data.frame(my_gene_col)
colnames(my_gene_col)[colnames(my_gene_col) == "my_gene_col"] <- "cluster"
my_gene_col$Cluster[my_gene_col$cluster == 1] <- "Cluster 1"
my_gene_col$Cluster[my_gene_col$cluster == 2] <- "Cluster 2"
my_gene_col$Cluster[my_gene_col$cluster == 3] <- "Cluster 3"
my_gene_col$Cluster[my_gene_col$cluster == 4] <- "Cluster 4"
my_gene_col$Cluster[my_gene_col$cluster == 5] <- "Cluster 5"
my_gene_col <- subset(my_gene_col, select = Cluster)
my_sample_col <- clindata[,c(8,7,2,33)]

DMPPredbetasHM <- pheatmap(DMPPredbetas_norm, annotation_row = my_gene_col, annotation_col = my_sample_col, cutree_cols = 3, fontsize_row = 4, fontsize_col = 5, 
                           cutree_rows = 5, show_rownames = F)

save_pheatmap_png <- function(x, filename, width=2400, height=2000, res = 300) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

save_pheatmap_png(DMPPredbetasHM, "DMPPredbetasHM.png")
