setwd("~/Documents/Code")

library("limma")
library("Biobase")
library(mixOmics)
library(rgl)
library(tidyverse)
library(ggpubr)

#Loading & tidying data----
load("~/Documents/Code/EXP.Rdata")

#Transpose transcriptomic data so individuals are in rows, genes in columns
Expdata_t <- (t(data))
rm(data)

#Remove 2 rows -9; Patient IDs #-----A and ----A
Expdata_t <- Expdata_t[-c(642,728),]

#Load clinical data with time from transcriptomic sample to breast stage 2 in years (TrT2)
ClinData_B2_Quartiles <- read.csv("~/Documents/Code/B2quartilesT.csv")

#Remove column X.1
ClinData_B2_Quartiles <- ClinData_B2_Quartiles[,-1]
rownames(ClinData_B2_Quartiles)<-paste(ClinData_B2_Quartiles[,1])
ClinData_B2_Quartiles <- ClinData_B2_Quartiles[,-1]

#Subset of transcriptomic data to include only those with a value for TrB2
Expdata_DtB2 <- Expdata_t[match(rownames(ClinData_B2_Quartiles),rownames(Expdata_t)),]
rm(Expdata_t)

#Check both files in same order
rownames(Expdata_DtB2)[c(10,65,275,361)]
rownames(ClinData_B2_Quartiles)[c(10,65,275,361)]

#Full B2 cohort before removal of outliers - PCA & PLSDA ----
#Quartile Groups for B2
B2Qgroups<- ClinData_B2_Quartiles[,c(2283,2292,2294)]
B2Qgroups$ClQGp <-B2Qgroups$DtoPQuartile
B2Qgroups$ClQGp[B2Qgroups$ClQGp=="2nd Quartile"] <- "2nd-4th Quartiles"
B2Qgroups$ClQGp[B2Qgroups$ClQGp=="3rd Quartile"] <- "2nd-4th Quartiles"
B2Qgroups$ClQGp[B2Qgroups$ClQGp=="4th Quartile"] <- "2nd-4th Quartiles"

B2Qgroups$FuQGp <-B2Qgroups$DtoPQuartile
B2Qgroups$FuQGp[B2Qgroups$FuQGp=="1st Quartile"] <- "1st-3rd Quartiles"
B2Qgroups$FuQGp[B2Qgroups$FuQGp=="2nd Quartile"] <- "1st-3rd Quartiles"
B2Qgroups$FuQGp[B2Qgroups$FuQGp=="3rd Quartile"] <- "1st-3rd Quartiles"

#Create object where B2 quartile is a factor
Y <- as.factor(B2Qgroups$DtoPQuartile)

#PCA all probes, B2
B2.pca <- pca(Expdata_DtB2) 
plotIndiv(B2.pca, group=Y, title="PCA, transcriptome, all probes", legend=TRUE, 
          legend.title="Time to B2", size.xlabel=rel(1.5), size.ylabel = rel(1.5), size.axis=rel(1.5), 
          size.legend.title=rel(1.5))

#B2 cohort (outliers removed) - PCA & PLSDA----
#Removing outliers identified on original PCA (all probes)
#Identify positions of outliers in clin data file
which(rownames(ClinData_B2_Quartiles) == "----A")
which(rownames(ClinData_B2_Quartiles) == "----A")
which(rownames(ClinData_B2_Quartiles) == "----A")
which(rownames(ClinData_B2_Quartiles) == "----A")
#etc.

#Create file with outliers removed
ClinData_B2_Quartiles_or <- ClinData_B2_Quartiles[-c(220, 267, 310, 309, 303, 192, 232, 178, 136, 208, 258, 343, 261, 285, 274, 129, 337, 177, 10, 211, 2, 157,266,315,320),]

#Recalculation of quartiles now that outliers removed
#Calculate quartiles of TrB2 & put in new object
QuartilesTrB2_or <- quantile(ClinData_B2_Quartiles_or$TrB2)

#Populate Quartile column with 1st-4th quartile (negative values included, outliers removed)
ClinData_B2_Quartiles_or$DtoPQuartile <- ifelse(ClinData_B2_Quartiles_or$TrB2<QuartilesTrB2_or[2],"1st Quartile", ClinData_B2_Quartiles_or$DtoPQuartile)
ClinData_B2_Quartiles_or$DtoPQuartile <- ifelse(ClinData_B2_Quartiles_or$TrB2>=QuartilesTrB2_or[2] & ClinData_B2_Quartiles_or$DtoPQuartile <QuartilesTrB2_or[3], "2nd Quartile", ClinData_B2_Quartiles_or$DtoPQuartile)
ClinData_B2_Quartiles_or$DtoPQuartile <- ifelse(ClinData_B2_Quartiles_or$TrB2>=QuartilesTrB2_or[3] & ClinData_B2_Quartiles_or$DtoPQuartile <QuartilesTrB2_or[4], "3rd Quartile", ClinData_B2_Quartiles_or$DtoPQuartile)
ClinData_B2_Quartiles_or$DtoPQuartile <- ifelse(ClinData_B2_Quartiles_or$TrB2>=QuartilesTrB2_or[4], "4th Quartile", ClinData_B2_Quartiles_or$DtoPQuartile)

#Subset of transcriptomic data to include only those participants in ClinData_B2_Quartiles_or
Expdata_DtB2_or <- Expdata_DtB2[match(rownames(ClinData_B2_Quartiles_or),rownames(Expdata_DtB2)),]

#Check both files in same order
rownames(Expdata_DtB2_or)[c(10,65,275,339)]
rownames(ClinData_B2_Quartiles_or)[c(10,65,275,339)]

#New PCA including all probes but with 25 outliers removed
#Quartile Groups for B2
B2Qgroups_or<- ClinData_B2_Quartiles_or[,c(2283,2292,2294)]
B2Qgroups_or$ClQGp <-B2Qgroups_or$DtoPQuartile
B2Qgroups_or$ClQGp[B2Qgroups_or$ClQGp=="2nd Quartile"] <- "2nd-4th Quartiles"
B2Qgroups_or$ClQGp[B2Qgroups_or$ClQGp=="3rd Quartile"] <- "2nd-4th Quartiles"
B2Qgroups_or$ClQGp[B2Qgroups_or$ClQGp=="4th Quartile"] <- "2nd-4th Quartiles"

B2Qgroups_or$FuQGp <-B2Qgroups_or$DtoPQuartile
B2Qgroups_or$FuQGp[B2Qgroups_or$FuQGp=="1st Quartile"] <- "1st-3rd Quartiles"
B2Qgroups_or$FuQGp[B2Qgroups_or$FuQGp=="2nd Quartile"] <- "1st-3rd Quartiles"
B2Qgroups_or$FuQGp[B2Qgroups_or$FuQGp=="3rd Quartile"] <- "1st-3rd Quartiles"

write.csv(B2Qgroups_or,file="B2Qgroups_or.csv")

#Create object where B2 quartile is a factor
YB2_or <- as.factor(B2Qgroups_or$DtoPQuartile)
YB2_or.pca <- pca(Expdata_DtB2_or) 
plotIndiv(YB2_or.pca, group=YB2_or, title="PCA, transcriptome, all probes, outliers removed", ind.names = FALSE, 
          legend=TRUE, legend.title="Time to B2", size.xlabel=rel(1.5), size.ylabel = rel(1.5), size.axis=rel(1.5), 
          size.legend.title=rel(1.5))

#PLSDA, all probes, B2
plsda.dataYB2_or<- plsda(Expdata_DtB2_or, YB2_or, ncomp = 339)

#Calculate variance
ev <- explained_variance(plsda.dataYB2_or$X, plsda.dataYB2_or$variates$X, ncomp =339)
Comp1v<-round(ev[1],digits = 2)*100
Comp2v<-round(ev[2],digits = 2)*100

plsdadata_plot <- data.frame(plsda.dataYB2_or$variates$X[,1:2])

ggplot(plsdadata_plot,aes(comp1,comp2,colour=YB2_or, shape=YB2_or))+
  geom_point(size=5, alpha=0.8)+
  coord_fixed(ratio=1)+
  theme_classic(base_size=30)+
  labs(title= "h) PLSDA, all probes", colour="Time to B2", shape="Time to B2")+
  xlab(label= paste("Component 1 (expl. variance= ", Comp1v, "%)",sep=""))+
  ylab(label= paste("Component 2 (expl. variance= ", Comp2v, "%)",sep=""))+
  theme(legend.position= "bottom")+
  scale_color_manual(values = c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728"), labels=c("Q1", "Q2", "Q3", "Q4"))+
  scale_shape_manual(values=c(1, 2, 3, 4), labels=c("Q1", "Q2", "Q3", "Q4"))+
  stat_stars(aes(colour=YB2_or), alpha=0.5)+
  stat_ellipse(aes(colour=YB2_or), level=0.95, linewidth = 1.5, type="t", alpha=0.5)
ggsave("SQhPLSDAall.png", width = 10, height = 10, dpi=320)

#B2 cohort (outliers removed) - DGE----
#Expression set - closest to B2 vs rest (outliers removed)
Eset_expdata_TrB2_1st_or <- ExpressionSet(t(Expdata_DtB2_or))
phenoData(Eset_expdata_TrB2_1st_or) <- AnnotatedDataFrame(B2Qgroups_or)
design_closesttB2_or <- model.matrix(~Eset_expdata_TrB2_1st_or@phenoData$ClQGp)
fit_closesttB2_or <- lmFit(Eset_expdata_TrB2_1st_or,design_closesttB2_or)
fit_closesttB2_or <- eBayes(fit_closesttB2_or)
list_closesttB2_or <- topTable(fit_closesttB2_or, number = 500, adjust.method = "BH")
write.csv(list_closesttB2_or,file="Top500_closesttB2_or.csv")

#Expression set - furthest from B2 vs rest (outliers removed)
Eset_expdata_TrB2_4th_or <- ExpressionSet(t(Expdata_DtB2_or))
phenoData(Eset_expdata_TrB2_4th_or) <- AnnotatedDataFrame(B2Qgroups_or)
design_furthestfB2_or <- model.matrix(~Eset_expdata_TrB2_4th_or@phenoData$FuQGp)
fit_furthestfB2_or  <- lmFit(Eset_expdata_TrB2_4th_or,design_furthestfB2_or)
fit_furthestfB2_or  <- eBayes(fit_furthestfB2_or)
list_furthestfB2_or  <- topTable(fit_furthestfB2_or, number = 500, adjust.method = "BH")
write.csv(list_furthestfB2_or ,file="Top500_furthestfB2_or.csv")

#PCA & PLSDA Analysis using top 500 differentially expressed probes for each group (outliers removed)
#Create vector of probe names for each group
Probes_ClB2_or <- rownames(list_closesttB2_or)
Probes_FuB2_or <- rownames(list_furthestfB2_or)

#Create subsets of transcription data: top 500 differentially expressed probes for each group (outliers removed)
Expdata_ClB2_or <- Expdata_DtB2_or[,match(Probes_ClB2_or,colnames(Expdata_DtB2_or))]
Expdata_FuB2_or <- Expdata_DtB2_or[,match(Probes_FuB2_or,colnames(Expdata_DtB2_or))]

write.csv(Expdata_ClB2_or,file="Expdata_ClB2_or.csv")
write.csv(Expdata_FuB2_or,file="Expdata_FuB2_or.csv")

#PLSDA, top 500 DGE, closest to B2 (outliers removed)
plsda.dataYClB2_or<- plsda(Expdata_ClB2_or, YB2_or, ncomp = 339)

#Calculate variance
ev <- explained_variance(plsda.dataYClB2_or$X, plsda.dataYClB2_or$variates$X, ncomp =339)
Comp1v<-round(ev[1],digits = 2)*100
Comp2v<-round(ev[2],digits = 2)*100

plsdadata_plot <- data.frame(plsda.dataYClB2_or$variates$X[,1:2])

ggplot(plsdadata_plot,aes(comp1,comp2,colour=YB2_or, shape=YB2_or))+
  geom_point(size=5, alpha=0.8)+
  coord_fixed(ratio=1)+
  theme_classic(base_size=30)+
  labs(title= "i) PLSDA, top 500, closest", colour="Time to B2", shape="Time to B2")+
  xlab(label= paste("Component 1 (expl. variance= ", Comp1v, "%)",sep=""))+
  ylab(label= paste("Component 2 (expl. variance= ", Comp2v, "%)",sep=""))+
  theme(legend.position= "bottom")+
  scale_color_manual(values = c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728"), labels=c("Q1", "Q2", "Q3", "Q4"))+
  scale_shape_manual(values=c(1, 2, 3, 4), labels=c("Q1", "Q2", "Q3", "Q4"))+
  stat_stars(aes(colour=YB2_or), alpha=0.5)+
  stat_ellipse(aes(colour=YB2_or), level=0.95, linewidth = 1.5, type="t", alpha=0.5)
ggsave("SQiPLSDA500cl.png", width = 10, height = 10, dpi=320)

#PCA top 500 DGE, closest to B2 (outliers removed)
YClB2_or.pca <- pca(Expdata_ClB2_or, ncomp = 339) 

#Calculate variance
eigs <- YClB2_or.pca$sdev^2
PCA1_var <- eigs[1] / sum(eigs)
PCA1_var<-round(PCA1_var,digits = 2)*100
PCA2_var <- eigs[2] / sum(eigs)
PCA2_var<-round(PCA2_var,digits = 2)*100

pcadata_plot <- data.frame(YClB2_or.pca$variates$X[,1:2])

ggplot(pcadata_plot,aes(PC1,PC2,colour=YB2_or, shape=YB2_or))+
  geom_point(size=5, alpha=0.8)+
  coord_fixed(ratio=1)+
  theme_classic(base_size=30)+
  labs(title= "f) PCA, top 500, closest", colour="Time to B2", shape="Time to B2")+
  xlab(label= paste("PC1 (expl. variance= ", PCA1_var, "%)",sep=""))+
  ylab(label= paste("PC2 (expl. variance= ", PCA2_var, "%)",sep=""))+
  theme(legend.position = "bottom")+
  scale_color_manual(values = c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728"), labels=c("Q1", "Q2", "Q3", "Q4"))+
  scale_shape_manual(values=c(1, 2, 3, 4), labels=c("Q1", "Q2", "Q3", "Q4"))+
  stat_stars(aes(colour=YB2_or), alpha=0.5)+
  stat_ellipse((aes(colour=YB2_or)), level=0.95, linewidth = 1.5, type="t", alpha=0.5)
ggsave("SQfPCA500cl.png", width = 10, height = 10, dpi=320)

#PLSDA top 500 DGE, furthest from B2 (outliers removed)
plsda.dataYFuB2_or<- plsda(Expdata_FuB2_or, YB2_or, ncomp = 339)

#Calculate variance
ev <- explained_variance(plsda.dataYFuB2_or$X, plsda.dataYFuB2_or$variates$X, ncomp =339)
Comp1v<-round(ev[1],digits = 2)*100
Comp2v<-round(ev[2],digits = 2)*100

plsdadata_plot <- data.frame(plsda.dataYFuB2_or$variates$X[,1:2])

ggplot(plsdadata_plot,aes(comp1,comp2,colour=YB2_or, shape=YB2_or))+
  geom_point(size=5, alpha=0.8)+
  coord_fixed(ratio=1)+
  theme_classic(base_size=30)+
  labs(title= "j) PLSDA, top 500, furthest", colour="Time to B2", shape="Time to B2")+
  xlab(label= paste("Component 1 (expl. variance= ", Comp1v, "%)",sep=""))+
  ylab(label= paste("Component 2 (expl. variance= ", Comp2v, "%)",sep=""))+
  theme(legend.position= "bottom")+
  scale_color_manual(values = c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728"), labels=c("Q1", "Q2", "Q3", "Q4"))+
  scale_shape_manual(values=c(1, 2, 3, 4), labels=c("Q1", "Q2", "Q3", "Q4"))+
  stat_stars(aes(colour=YB2_or), alpha=0.5)+
  stat_ellipse((aes(colour=YB2_or)), level=0.95, linewidth = 1.5, type="t", alpha=0.5)
ggsave("SQjPLSDA500fu.png", width = 10, height = 10, dpi=320)

#PCA top 500 DGE, furthest from B2 (outliers removed)
YFuB2_or.pca <- pca(Expdata_FuB2_or, ncomp = 339) 

#Calculate variance
eigs <- YFuB2_or.pca$sdev^2
PCA1_var <- eigs[1] / sum(eigs)
PCA1_var<-round(PCA1_var,digits = 2)*100

PCA2_var <- eigs[2] / sum(eigs)
PCA2_var<-round(PCA2_var,digits = 2)*100

pcadata_plot <- data.frame(YFuB2_or.pca$variates$X[,1:2])

ggplot(pcadata_plot,aes(PC1,PC2,colour=YB2_or, shape=YB2_or))+
  geom_point(size=5, alpha=0.8)+
  coord_fixed(ratio=1)+
  theme_classic(base_size=30)+
  labs(title= "g) PCA, top 500,furthest", colour="Time to B2", shape="Time to B2")+
  xlab(label= paste("PC1 (expl. variance= ", PCA1_var, "%)",sep=""))+
  ylab(label= paste("PC2 (expl. variance= ", PCA2_var, "%)",sep=""))+
  theme(legend.position = "bottom")+
  scale_color_manual(values = c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728"), labels=c("Q1", "Q2", "Q3", "Q4"))+
  scale_shape_manual(values=c(1, 2, 3, 4), labels=c("Q1", "Q2", "Q3", "Q4"))+
  stat_stars(aes(colour=YB2_or), alpha=0.5)+
  stat_ellipse((aes(colour=YB2_or)), level=0.95, linewidth = 1.5, type="t", alpha=0.5)
ggsave("SQgPCA500fu.png", width = 10, height = 10, dpi=320)

