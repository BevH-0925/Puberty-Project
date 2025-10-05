library(readr)
library(edgeR)
library(stringr)
library(EnsDb.Hsapiens.v79)
library(mixOmics)
library(limma)
library(tidyverse)
library(ggpubr)

#Load Data----
setwd('~/Documents/RNASeqBH1-BH32/EdgeR')
data <- read.table('~/Documents/RNASeqBH1-BH32/EdgeR/Data/bev_counts.txt', sep = "\t", header = TRUE)

#Add in gene names
rownames(data) <- data$Geneid 
gene_names <- ensembldb::select(EnsDb.Hsapiens.v79, keys= row.names(data), keytype = "GENEID", columns = c("GENENAME","GENEID"))# changes the ENSEMBL names to gene names
gene_names$GENENAME <- make.unique(gene_names$GENENAME,sep = "~")

data <- data[row.names(data) %in% gene_names$GENEID,]#remove any genes that don't match with identifier 

#Check if data$Geneid is in the same order as gene_names$GENEID
is_same_order <- all(identical(data$Geneid, gene_names$GENEID))
#Print the result
print(is_same_order)

#Merge the data frames based on Geneid column
merged_data <- merge(data, gene_names, by.x = "Geneid", by.y = "GENEID", all.x = TRUE)
#Move GENENAME column to the beginning
merged_data <- cbind(merged_data[ncol(merged_data)], merged_data[-ncol(merged_data)])

data1 <-merged_data

drop <- c("Chr", "Start", "End", "Strand", "Length", "Geneid") # removing stuff you don't need
data1 <- data1[, !(names(data1) %in% drop)]
rownames(data1) <- data1$GENENAME

#Changes sample names to shorten them
sample_names <- strsplit(colnames(data1),"\\.")
sample_names <- unlist(lapply(sample_names[-1], "[[", 6))# select the element you want ignoring the first columns

#set col names of data
colnames(data1) <- c("GENENAME", sample_names)

sample_info <- read_csv("~/Documents/RNASeqBH1-BH32/Clinical Data Analysis/Data/BH1-32.csv")

#Reorder count data cols to match the sample_info data or the way around
data1 <- data1[,c(1, na.omit(match(sample_info$TranscriptomeID,colnames(data1))))]

#Tidy up environment
rm(data)
rm(merged_data)
rm(is_same_order)
rm(drop)
rm(sample_names)

sample_info$Ethnicity <- as.factor(sample_info$Ethnicity)
sample_info$Group <- as.factor(sample_info$Group)
sample_info$Sex <- as.factor(sample_info$Sex)
sample_info$AgeGroups<- ifelse(sample_info$Age_yrs<12, "2to11years", "14to18years")
sample_info$AgeGroups <- as.factor(sample_info$AgeGroups)

#Create comparison groups----
group <- sample_info$Group
ethnicity <- sample_info$Ethnicity
sex <- sample_info$Sex
agegp <- sample_info$AgeGroups
design <- model.matrix(~0+group+ethnicity+sex+agegp)
my_contrast <- makeContrasts(groupPubertal - groupPrepubertal,levels= design)

#Filter and normalise data----
edgeRformat_data <- DGEList(counts=data1, group=group)
keep <- filterByExpr(edgeRformat_data) #ID low expression genes
edgeRformat_data <- edgeRformat_data[keep,,keep.lib.sizes=F] #remove low expression genes
edgeRformat_data <- normLibSizes(edgeRformat_data) #normalises for library sizes
fit <- glmQLFit(edgeRformat_data, design, robust=TRUE)
plotQLDisp(fit)

#Run DGE (glm method)----
qlf<- glmQLFTest(fit,contrast=my_contrast)
summary(decideTests(qlf))
plotMD(qlf)
dge_output <- topTags(qlf, n=500, sort.by = "PValue")
write.csv(dge_output, row.names=TRUE, "~/Documents/RNASeqBH1-BH32/EdgeR/DEGs500_BH1-32T_glmsea.csv")
dge_output_all <- topTags(qlf, n=19585, sort.by = "PValue")
write.csv(dge_output_all,row.names=TRUE, "~/Documents/RNASeqBH1-BH32/EdgeR/DEGsall_BH1-32T_glmsea.csv")

#Top 500 DEGs----
logcounts <- cpm(edgeRformat_data,log = T)
lcpm.df <- as.data.frame(logcounts)
write.csv(lcpm.df, "~/Documents/RNASeqBH1-BH32/EdgeR/BH1-32T_lcpm.csv")

# Check distributions of samples using boxplots
boxplot(logcounts, xlab="Sample Name", ylab="Log2 counts per million",las=2) #las rotates sample name
# Add a blue horizontal line that corresponds to the median logCPM
abline(h=median(logcounts),col="blue")
title("Boxplots of logCPMs")

top_500 <- dge_output[["table"]]
top500_lcpm <- logcounts[na.omit(match(rownames(top_500),rownames(logcounts))),]
top500_lcpm.df <- as.data.frame(top500_lcpm)
write.csv(top500_lcpm.df, "~/Documents/RNASeqBH1-BH32/EdgeR/BH1-32T_top500_lcpm_glmsea.csv")

#PCA all MF ##----
pcadata_all<-pca(t(logcounts),ncomp=32)
Y <- as.factor(group)

#Calculate variance
eigs <- pcadata_all$sdev^2
PCA1_all_var <- eigs[1] / sum(eigs)
PCA1_all_var<-round(PCA1_all_var,digits = 2)*100

PCA2_all_var <- eigs[2] / sum(eigs)
PCA2_all_var<-round(PCA2_all_var,digits = 2)*100

pcadata_plot <- data.frame(pcadata_all$variates$X[,1:2])

#Note default type t distribution used for ellipse as better for small samples. 95% confidence ellipse.
ggplot(pcadata_plot,aes(PC1,PC2,colour=Y, shape=Y))+
  geom_point(size=7)+
  coord_fixed(ratio=1)+
  theme_classic(base_size=30)+
  labs(title= "a) PCA, all genes, MF", colour="Pubertal Status", shape="Pubertal Status")+
  xlab(label= paste("PC1 (expl. variance= ", PCA1_all_var, "%)",sep=""))+
  ylab(label= paste("PC2 (expl. variance= ", PCA2_all_var, "%)",sep=""))+
  theme(legend.position = "bottom")+
  scale_color_manual(values = c("#1F77B4", "#FF7F0E"))+
  stat_stars(aes(colour=Y))+
  stat_ellipse((aes(colour=Y)), level=0.95, linewidth = 1.5, type="t")
ggsave("SQaPCAallgenesMFpuberty.png", width = 10, height = 10, dpi=320)
  
Z <- as.factor(sample_info$Sex)
ggplot(pcadata_plot,aes(PC1,PC2,colour=Z, shape=Z))+
  geom_point(size=7)+
  coord_fixed(ratio=1)+
  theme_classic(base_size=30)+
  labs(title= "b) PCA, all genes, MF", colour="Sex", shape="Sex")+
  xlab(label= paste("PC1 (expl. variance= ", PCA1_all_var, "%)",sep=""))+
  ylab(label= paste("PC2 (expl. variance= ", PCA2_all_var, "%)",sep=""))+
  theme(legend.position = "bottom")+
  scale_color_manual(values = c("#2CA02C", "#D62728"))+
  stat_stars(aes(colour=Z))+
  stat_ellipse((aes(colour=Z)), level=0.95, linewidth = 1.5, type="t")
ggsave("SQbPCAallgenesMFsex.png", width = 10, height = 10, dpi=320)

A <- as.factor(sample_info$Ethnicity)
ggplot(pcadata_plot,aes(PC1,PC2,colour=A, shape=A))+
  geom_point(size=7)+
  coord_fixed(ratio=1)+
  theme_classic(base_size=30)+
  labs(title= "c) PCA, all genes, MF", colour="Ethnic Group", shape="Ethnic Group")+
  xlab(label= paste("PC1 (expl. variance= ", PCA1_all_var, "%)",sep=""))+
  ylab(label= paste("PC2 (expl. variance= ", PCA2_all_var, "%)",sep=""))+
  theme(legend.position = "bottom")+
  scale_color_manual(values = c("#9467BD", "#8C564B"))+
  stat_stars(aes(colour=A))+
  stat_ellipse((aes(colour=A)), level=0.95, linewidth = 1.5, type="t")
ggsave("SQcPCAallgenesMFethnicity.png", width = 10, height = 10, dpi=320)

B <- as.factor(sample_info$AgeGroups)
ggplot(pcadata_plot,aes(PC1,PC2,colour=B, shape=B))+
  geom_point(size=7)+
  coord_fixed(ratio=1)+
  theme_classic(base_size=30)+
  labs(title= "d) PCA, all genes, MF", colour="Age Group", shape="Age Group")+
  xlab(label= paste("PC1 (expl. variance= ", PCA1_all_var, "%)",sep=""))+
  ylab(label= paste("PC2 (expl. variance= ", PCA2_all_var, "%)",sep=""))+
  theme(legend.position = "bottom")+
  scale_color_manual(values = c("#E377C2", "#7F7F7F"))+
  stat_stars(aes(colour=B))+
  stat_ellipse((aes(colour=B)), level=0.95, linewidth = 1.5, type="t")
ggsave("SQdPCAallgenesMFagegpNTA.png", width = 10, height = 10, dpi=320)

#Explained variance
plot(pcadata_all)
plotLoadings(pcadata_all)

#PLSDA all MF----
plsda_all<- plsda(t(logcounts), Y, ncomp = 32)

#Calculate variance
ev <- explained_variance(plsda_all$X, plsda_all$variates$X, ncomp =32)
Comp1v<-round(ev[1],digits = 2)*100
Comp2v<-round(ev[2],digits = 2)*100

plsdadata_plot <- data.frame(plsda_all$variates$X[,1:2])

ggplot(plsdadata_plot,aes(comp1,comp2,colour=Y, shape=Y))+
  geom_point(size=7)+
  coord_fixed(ratio = 1)+
  theme_classic(base_size=30)+
  labs(title= "e) PLSDA, all genes, MF", colour="Pubertal Status", shape="Pubertal Status")+
  xlab(label= paste("Component 1 (expl. variance= ", Comp1v, "%)",sep=""))+
  ylab(label= paste("Component 2 (expl. variance= ", Comp2v, "%)",sep=""))+
  theme(legend.position= "bottom")+
  scale_color_manual(values = c("#1F77B4", "#FF7F0E"))+
  stat_stars(aes(colour=Y))+
  stat_ellipse((aes(colour=Y)), level=0.95, linewidth = 1.5, type="t")
ggsave("SQePLSDAallgenesMFpubertyNTA.png", width = 10, height = 10, dpi=320)

plsda_all<- plsda(t(logcounts), Z, ncomp = 32)
ev <- explained_variance(plsda_all$X, plsda_all$variates$X, ncomp =32)
Comp1v<-round(ev[1],digits = 2)*100
Comp2v<-round(ev[2],digits = 2)*100
plsdadata_plot <- data.frame(plsda_all$variates$X[,1:2])

ggplot(plsdadata_plot,aes(comp1,comp2,colour=Z, shape=Z))+
  geom_point(size=7)+
  coord_fixed(ratio = 1)+
  theme_classic(base_size=30)+
  labs(title= "f) PLSDA, all genes, MF", colour="Sex", shape="Sex")+
  xlab(label= paste("Component 1 (expl. variance= ", Comp1v, "%)",sep=""))+
  ylab(label= paste("Component 2 (expl. variance= ", Comp2v, "%)",sep=""))+
  theme(legend.position= "bottom")+
  scale_color_manual(values = c("#2CA02C", "#D62728"))+
  stat_stars(aes(colour=Z))+
  stat_ellipse((aes(colour=Z)), level=0.95, linewidth = 1.5, type="t")
ggsave("SQfPLSDAallgenesMFsexNTA.png", width = 10, height = 10, dpi=320)

plsda_all<- plsda(t(logcounts), A, ncomp = 32)
ev <- explained_variance(plsda_all$X, plsda_all$variates$X, ncomp =32)
Comp1v<-round(ev[1],digits = 2)*100
Comp2v<-round(ev[2],digits = 2)*100
plsdadata_plot <- data.frame(plsda_all$variates$X[,1:2])

ggplot(plsdadata_plot,aes(comp1,comp2,colour=A, shape=A))+
  geom_point(size=7)+
  coord_fixed(ratio = 1)+
  theme_classic(base_size=30)+
  labs(title= "g) PLSDA, all genes, MF", colour="Ethnic Group", shape="Ethnic Group")+
  xlab(label= paste("Component 1 (expl. variance= ", Comp1v, "%)",sep=""))+
  ylab(label= paste("Component 2 (expl. variance= ", Comp2v, "%)",sep=""))+
  theme(legend.position= "bottom")+
  scale_color_manual(values = c("#9467BD", "#8C564B"))+
  stat_stars(aes(colour=A))+
  stat_ellipse((aes(colour=A)), level=0.95, linewidth = 1.5, type="t")
ggsave("SQgPLSDAallgenesMFethnicity.png", width = 10, height = 10, dpi=320)


plsda_all<- plsda(t(logcounts), B, ncomp = 32)
ev <- explained_variance(plsda_all$X, plsda_all$variates$X, ncomp =32)
Comp1v<-round(ev[1],digits = 2)*100
Comp2v<-round(ev[2],digits = 2)*100
plsdadata_plot <- data.frame(plsda_all$variates$X[,1:2])

ggplot(plsdadata_plot,aes(comp1,comp2,colour=B, shape=B))+
  geom_point(size=7)+
  coord_fixed(ratio = 1)+
  theme_classic(base_size=30)+
  labs(title= "h) PLSDA, all genes, MF", colour="Age Group", shape="Age Group")+
  xlab(label= paste("Component 1 (expl. variance= ", Comp1v, "%)",sep=""))+
  ylab(label= paste("Component 2 (expl. variance= ", Comp2v, "%)",sep=""))+
  theme(legend.position= "bottom")+
  scale_color_manual(values = c("#E377C2", "#7F7F7F"))+
  stat_stars(aes(colour=B))+
  stat_ellipse((aes(colour=B)), level=0.95, linewidth = 1.5, type="t")
ggsave("SQhPLSDAallgenesMFagegpNTA.png", width = 10, height = 10, dpi=320)

#PCA top 500 - contains normalised counts for all the data MF----
pcadata_500<- pca(t(top500_lcpm), ncomp=32)
Y <- as.factor(group)

#Calculate variance
eigs <- pcadata_500$sdev^2
PCA1_500_var <- eigs[1] / sum(eigs)
PCA1_500_var<-round(PCA1_500_var,digits = 2)*100

eigs <- pcadata_500$sdev^2
PCA2_500_var <- eigs[2] / sum(eigs)
PCA2_500_var<-round(PCA2_500_var,digits = 2)*100

pcadata_plot <- data.frame(pcadata_500$variates$X[,1:2])

ggplot(pcadata_plot,aes(PC1,PC2,colour=Y, shape=Y))+
  geom_point(size=7)+
  coord_fixed(ratio=1)+
  theme_classic(base_size=30)+
  labs(title= "b) PCA, top 500, glm, MF", colour="Pubertal Status", shape="Pubertal Status")+
  xlab(label= paste("PC1 (expl. variance= ", PCA1_500_var, "%)",sep=""))+
  ylab(label= paste("PC2 (expl. variance= ", PCA2_500_var, "%)",sep=""))+
  theme(legend.position = "bottom")+
  scale_color_manual(values = c("#1F77B4", "#FF7F0E"))+
  stat_stars(aes(colour=Y))+
  stat_ellipse((aes(colour=Y)), level=0.95, linewidth = 1.5, type="t")
ggsave("SQbPCAtop500MFpubertyglm.png", width = 10, height = 10, dpi=320)

#Explained variance
plot(pcadata_500)
plotLoadings(pcadata_500)

#PLSDA top 500 MF----
plsda_500<- plsda(t(top500_lcpm), Y, ncomp = 32)

#Calculate variance
ev <- explained_variance(plsda_500$X, plsda_500$variates$X, ncomp =32)
Comp1v<-round(ev[1],digits = 2)*100
Comp2v<-round(ev[2],digits = 2)*100

plsdadata_plot <- data.frame(plsda_500$variates$X[,1:2])

ggplot(plsdadata_plot,aes(comp1,comp2,colour=Y, shape=Y))+
  geom_point(size=7)+
  coord_fixed(ratio = 1)+
  theme_classic(base_size=30)+
  labs(title= "d) PLSDA, top 500, glm, MF", colour="Pubertal Status", shape="Pubertal Status")+
  xlab(label= paste("Component 1 (expl. variance= ", Comp1v, "%)",sep=""))+
  ylab(label= paste("Component 2 (expl. variance= ", Comp2v, "%)",sep=""))+
  theme(legend.position= "bottom")+
  scale_color_manual(values = c("#1F77B4", "#FF7F0E"))+
  stat_stars(aes(colour=Y))+
  stat_ellipse((aes(colour=Y)), level=0.95, linewidth = 1.5, type="t")
ggsave("SQdPLSDAtop500MFpubertyglm.png", width = 10, height = 10, dpi=320)

#Classic edgeR Pipeline MF----
group <- sample_info$Group
edgeRformat_data <- DGEList(counts=data1, group=group)
keep <- filterByExpr(edgeRformat_data) #ID low expression genes
edgeRformat_data <- edgeRformat_data[keep,,keep.lib.sizes=F]#remove low expression genes
edgeRformat_data <- normLibSizes(edgeRformat_data) #normalises for library sizes
edgeRformat_data <- estimateDisp(edgeRformat_data, robust=TRUE)
plotBCV(edgeRformat_data)#plot of biological CV
et <- exactTest(edgeRformat_data)
dge_output_classic <- topTags(et, n=500, sort.by = "PValue")
write.csv(dge_output_classic, row.names=TRUE, "~/Documents/Files from PEC drive/RNASeqBH1-BH32/EdgeR/DEGs500_BH1-32T_classic.csv")
dge_output_classic_all <- topTags(et, n=19585, sort.by = "PValue")
write.csv(dge_output_classic_all,row.names=TRUE, "~/Documents/Files from PEC drive/RNASeqBH1-BH32/EdgeR/DEGsall_BH1-32T_classic.csv")

logcounts <- cpm(edgeRformat_data,log = T)

# Check distributions of samples using boxplots
boxplot(logcounts, xlab="Sample Name", ylab="Log2 counts per million",las=2) #las rotates sample name
# Add a blue horizontal line that corresponds to the median logCPM
abline(h=median(logcounts),col="blue")
title("Boxplots of logCPMs")

top_500_classic <- dge_output_classic[["table"]]
top500_lcpm_classic <- logcounts[na.omit(match(rownames(top_500_classic),rownames(logcounts))),]
top500_lcpm.df_classic <- as.data.frame(top500_lcpm_classic)
write.csv(top500_lcpm.df_classic, "~/Documents/RNASeqBH1-BH32/EdgeR/BH1-32T_top500_lcpm_classic.csv")

#PCA classic MF----
pcadata_all<-pca(t(logcounts),ncomp=32)
Y <- as.factor(group)

# PCA of top 500 - contains normalised counts for all the data
pcadata_500<- pca(t(top500_lcpm_classic), ncomp=32)

#Calculate variance
eigs <- pcadata_500$sdev^2
PCA1_500_var <- eigs[1] / sum(eigs)
PCA1_500_var<-round(PCA1_500_var,digits = 2)*100

eigs <- pcadata_500$sdev^2
PCA2_500_var <- eigs[2] / sum(eigs)
PCA2_500_var<-round(PCA2_500_var,digits = 2)*100

pcadata_plot <- data.frame(pcadata_500$variates$X[,1:2])

ggplot(pcadata_plot,aes(PC1,PC2,colour=Y, shape=Y))+
  geom_point(size=7)+
  coord_fixed(ratio=1)+
  theme_classic(base_size=30)+
  labs(title= "a) PCA, top 500, classic, MF", colour="Pubertal Status", shape="Pubertal Status")+
  xlab(label= paste("PC1 (expl. variance= ", PCA1_500_var, "%)",sep=""))+
  ylab(label= paste("PC2 (expl. variance= ", PCA2_500_var, "%)",sep=""))+
  theme(legend.position = "bottom")+
  scale_color_manual(values = c("#1F77B4", "#FF7F0E"))+
  stat_stars(aes(colour=Y))+
  stat_ellipse((aes(colour=Y)), level=0.95, linewidth = 1.5, type="t")
ggsave("SQaPCAtop500MFpubertyclassic.png", width = 10, height = 10, dpi=320)

#Explained variance
plot(pcadata_500)
plotLoadings(pcadata_500)

#PLSDA top 500 classic MF----
plsda_500<- plsda(t(top500_lcpm_classic), Y, ncomp = 32)

#Calculate variance
ev <- explained_variance(plsda_500$X, plsda_500$variates$X, ncomp =32)
Comp1v<-round(ev[1],digits = 2)*100
Comp2v<-round(ev[2],digits = 2)*100

plsdadata_plot <- data.frame(plsda_500$variates$X[,1:2])

ggplot(plsdadata_plot,aes(comp1,comp2,colour=Y, shape=Y))+
  geom_point(size=7)+
  coord_fixed(ratio=1)+
  theme_classic(base_size=30)+
  labs(title= "c) PLSDA, top 500, classic, MF", colour="Pubertal Status", shape="Pubertal Status")+
  xlab(label= paste("Component 1 (expl. variance= ", Comp1v, "%)",sep=""))+
  ylab(label= paste("Component 2 (expl. variance= ", Comp2v, "%)",sep=""))+
  theme(legend.position= "bottom")+
  scale_color_manual(values = c("#1F77B4", "#FF7F0E"))+
  stat_stars(aes(colour=Y))+
  stat_ellipse((aes(colour=Y)), level=0.95, linewidth = 1.5, type="t")
ggsave("SQcPLSDAtop500MFpubertyclassicNTA.png", width = 10, height = 10, dpi=320)

#Redo analysis for females only----
sample_info_females <- sample_info[sample_info$Sex=="Female",]
data1_females <- data1[,c(1, na.omit(match(sample_info_females$TranscriptomeID,colnames(data1))))]

#Create comparison groups - females
group <- sample_info_females$Group
ethnicity <- sample_info_females$Ethnicity
agegp <- sample_info_females$AgeGroups
design <- model.matrix(~0+group+ethnicity+agegp)
my_contrast <- makeContrasts(groupPubertal - groupPrepubertal,levels= design)

#Filter and normalise data F----
edgeRformat_data <- DGEList(counts=data1_females, group=group)
keep <- filterByExpr(edgeRformat_data) #ID low expression genes
edgeRformat_data <- edgeRformat_data[keep,,keep.lib.sizes=F] #remove low expression genes
edgeRformat_data <- normLibSizes(edgeRformat_data) #normalises for library sizes
fit <- glmQLFit(edgeRformat_data, design, robust=TRUE)
plotQLDisp(fit)

#Run DGE (glm method) F----
qlf<- glmQLFTest(fit,contrast=my_contrast)
summary(decideTests(qlf))
plotMD(qlf)
dge_output_F <- topTags(qlf, n=500, sort.by = "PValue")
write.csv(dge_output_F, row.names=TRUE, "~/Documents/RNASeqBH1-BH32/EdgeR/DEGs500_BH1-32T_F_glmea.csv")
dge_output_all_F <- topTags(qlf, n=19748, sort.by = "PValue")
write.csv(dge_output_all_F,row.names=TRUE, "~/Documents/RNASeqBH1-BH32/EdgeR/DEGsall_BH1-32T_F_glmea.csv")

#Top 500 DEGs F----
logcounts_F <- cpm(edgeRformat_data,log = T)
lcpm.df_F <- as.data.frame(logcounts_F)
write.csv(lcpm.df_F, "~/Documents/RNASeqBH1-BH32/EdgeR/BH1-32T_lcpm_F.csv")

# Check distributions of samples using boxplots
boxplot(logcounts_F, xlab="Sample Name", ylab="Log2 counts per million",las=2) #las rotates sample name
# Add a blue horizontal line that corresponds to the median logCPM
abline(h=median(logcounts_F),col="blue")
title("Boxplots of logCPMs")

top_500_F <- dge_output_F[["table"]]
top500_lcpm_F <- logcounts_F[na.omit(match(rownames(top_500_F),rownames(logcounts_F))),]
top500_lcpm_F.df <- as.data.frame(top500_lcpm_F)
write.csv(top500_lcpm_F.df, "~/Documents/RNASeqBH1-BH32/EdgeR/BH1-32T_top500_lcpm_F_glmea.csv")

#PCA all F glm----
pcadata_all_F<-pca(t(logcounts_F),ncomp=21)
Y <- as.factor(group)

#Calculate variance
eigs <- pcadata_all_F$sdev^2
PCA1_all_var <- eigs[1] / sum(eigs)
PCA1_all_var<-round(PCA1_all_var,digits = 2)*100
PCA2_all_var <- eigs[2] / sum(eigs)
PCA2_all_var<-round(PCA2_all_var,digits = 2)*100

pcadata_plot <- data.frame(pcadata_all_F$variates$X[,1:2])

ggplot(pcadata_plot,aes(PC1,PC2,colour=Y, shape=Y))+
  geom_point(size=7)+
  coord_fixed(ratio=1)+
  theme_classic(base_size=30)+
  labs(title= "a) PCA, all genes, F", colour="Pubertal Status", shape="Pubertal Status")+
  xlab(label= paste("PC1 (expl. variance= ", PCA1_all_var, "%)",sep=""))+
  ylab(label= paste("PC2 (expl. variance= ", PCA2_all_var, "%)",sep=""))+
  theme(legend.position = "bottom")+
  scale_color_manual(values = c("#1F77B4", "#FF7F0E"))+
  stat_stars(aes(colour=Y))+
  stat_ellipse((aes(colour=Y)), level=0.95, linewidth = 1.5, type="t")
ggsave("SQaPCAallgenesFpubertyNTA.png", width = 10, height = 10, dpi=320)

A <- as.factor(sample_info_females$Ethnicity)

ggplot(pcadata_plot,aes(PC1,PC2,colour=A, shape=A))+
  geom_point(size=7)+
  coord_fixed(ratio=1)+
  theme_classic(base_size=30)+
  labs(title= "b) PCA, all genes, F", colour="Ethnic Group", shape="Ethnic Group")+
  xlab(label= paste("PC1 (expl. variance= ", PCA1_all_var, "%)",sep=""))+
  ylab(label= paste("PC2 (expl. variance= ", PCA2_all_var, "%)",sep=""))+
  theme(legend.position = "bottom")+
  scale_color_manual(values = c("#9467BD", "#8C564B"))+
  stat_stars(aes(colour=A))+
  stat_ellipse((aes(colour=A)), level=0.95, linewidth = 1.5, type="t")
ggsave("SQbPCAallgenesFethnicity.png", width = 10, height = 10, dpi=320)

#Explained variance
plot(pcadata_all_F)
plotLoadings(pcadata_all_F)

#PLSDA all Fglm----
plsda_all_F<- plsda(t(logcounts_F), Y, ncomp = 21)

#Calculate variance
ev <- explained_variance(plsda_all_F$X, plsda_all_F$variates$X, ncomp =21)
Comp1v<-round(ev[1],digits = 2)*100
Comp2v<-round(ev[2],digits = 2)*100

plsdadata_plot <- data.frame(plsda_all_F$variates$X[,1:2])

ggplot(plsdadata_plot,aes(comp1,comp2,colour=Y, shape=Y))+
  geom_point(size=7)+
  coord_fixed(ratio = 1)+
  theme_classic(base_size=30)+
  labs(title= "c) PLSDA, all genes, F", colour="Pubertal Status", shape="Pubertal Status")+
  xlab(label= paste("Component 1 (expl. variance= ", Comp1v, "%)",sep=""))+
  ylab(label= paste("Component 2 (expl. variance= ", Comp2v, "%)",sep=""))+
  theme(legend.position= "bottom")+
  scale_color_manual(values = c("#1F77B4", "#FF7F0E"))+
  stat_stars(aes(colour=Y))+
  stat_ellipse((aes(colour=Y)), level=0.95, linewidth = 1.5, type="t")
ggsave("SQcPLSDAallgenesFpubertyNTA.png", width = 10, height = 10, dpi=320)

plsda_all_F<- plsda(t(logcounts_F), A, ncomp = 21)
ev <- explained_variance(plsda_all_F$X, plsda_all_F$variates$X, ncomp =21)
Comp1v<-round(ev[1],digits = 2)*100
Comp2v<-round(ev[2],digits = 2)*100
plsdadata_plot <- data.frame(plsda_all_F$variates$X[,1:2])

ggplot(plsdadata_plot,aes(comp1,comp2,colour=A, shape=A))+
  geom_point(size=7)+
  coord_fixed(ratio = 1)+
  theme_classic(base_size=30)+
  labs(title= "d) PLSDA, all genes, F", colour="Ethnic Group", shape="Ethnic Group")+
  xlab(label= paste("Component 1 (expl. variance= ", Comp1v, "%)",sep=""))+
  ylab(label= paste("Component 2 (expl. variance= ", Comp2v, "%)",sep=""))+
  theme(legend.position= "bottom")+
  scale_color_manual(values = c("#9467BD", "#8C564B"))+
  stat_stars(aes(colour=A))+
  stat_ellipse((aes(colour=A)), level=0.95, linewidth = 1.5, type="t")
ggsave("SQdPLSDAallgenesFethnicityNTA.png", width = 10, height = 10, dpi=320)

#PCA top 500 Fglm----
pcadata_500_F<- pca(t(top500_lcpm_F), ncomp=21)

#Calculate variance
eigs <- pcadata_500_F$sdev^2
PCA1_500_var <- eigs[1] / sum(eigs)
PCA1_500_var<-round(PCA1_500_var,digits = 2)*100
PCA2_500_var <- eigs[2] / sum(eigs)
PCA2_500_var<-round(PCA2_500_var,digits = 2)*100

pcadata_plot <- data.frame(pcadata_500_F$variates$X[,1:2])

ggplot(pcadata_plot,aes(PC1,PC2,colour=Y, shape=Y))+
  geom_point(size=7)+
  coord_fixed(ratio=1)+
  theme_classic(base_size=30)+
  labs(title= "f) PCA, top 500, glm, F", colour="Pubertal Status", shape="Pubertal Status")+
  xlab(label= paste("PC1 (expl. variance= ", PCA1_500_var, "%)",sep=""))+
  ylab(label= paste("PC2 (expl. variance= ", PCA2_500_var, "%)",sep=""))+
  theme(legend.position = "bottom")+
  scale_color_manual(values = c("#1F77B4", "#FF7F0E"))+
  stat_stars(aes(colour=Y))+
  stat_ellipse((aes(colour=Y)), level=0.95, linewidth = 1.5, type="t")
ggsave("SQfPCAtop500genesFpubertyglm.png", width = 10, height = 10, dpi=320)

#Explained variance
plot(pcadata_500_F)
plotLoadings(pcadata_500_F)

#PLSDA of top 500
plsda_500_F<- plsda(t(top500_lcpm_F), Y, ncomp = 21)
#Calculate variance
ev <- explained_variance(plsda_500_F$X, plsda_500_F$variates$X, ncomp =21)
Comp1v<-round(ev[1],digits = 2)*100
Comp2v<-round(ev[2],digits = 2)*100

plsdadata_plot <- data.frame(plsda_500_F$variates$X[,1:2])

ggplot(plsdadata_plot,aes(comp1,comp2,colour=Y, shape=Y))+
  geom_point(size=7)+
  coord_fixed(ratio = 1)+
  theme_classic(base_size=30)+
  labs(title= "h) PLSDA, top 500, glm, F", colour="Pubertal Status", shape="Pubertal Status")+
  xlab(label= paste("Component 1 (expl. variance= ", Comp1v, "%)",sep=""))+
  ylab(label= paste("Component 2 (expl. variance= ", Comp2v, "%)",sep=""))+
  theme(legend.position= "bottom")+
  scale_color_manual(values = c("#1F77B4", "#FF7F0E"))+
  stat_stars(aes(colour=Y))+
  stat_ellipse((aes(colour=Y)), level=0.95, linewidth = 1.5, type="t")
ggsave("SQhPLSDAtop500genesFpubertyglm.png", width = 10, height = 10, dpi=320)

#Classic edgeR Pipeline F----
group <- sample_info_females$Group
edgeRformat_data <- DGEList(counts=data1_females, group=group)
keep <- filterByExpr(edgeRformat_data) #ID low expression genes
edgeRformat_data <- edgeRformat_data[keep,,keep.lib.sizes=F]#remove low expression genes
edgeRformat_data <- normLibSizes(edgeRformat_data) #normalises for library sizes
edgeRformat_data <- estimateDisp(edgeRformat_data, robust=TRUE)
plotBCV(edgeRformat_data)#plot of biological CV
et <- exactTest(edgeRformat_data)
dge_output_classic_F <- topTags(et, n=500, sort.by = "PValue")
write.csv(dge_output_classic_F, row.names=TRUE, "~/Documents/RNASeqBH1-BH32/EdgeR/DEGs500_BH1-32T_classic_F.csv")
dge_output_classic_all_F <- topTags(et, n=19748, sort.by = "PValue")
write.csv(dge_output_classic_all_F,row.names=TRUE, "~/Documents/RNASeqBH1-BH32/EdgeR/DEGsall_BH1-32T_classic_F.csv")

logcounts_F<- cpm(edgeRformat_data,log = T)

# Check distributions of samples using boxplots
boxplot(logcounts_F, xlab="Sample Name", ylab="Log2 counts per million",las=2) #las rotates sample name
# Add a blue horizontal line that corresponds to the median logCPM
abline(h=median(logcounts_F),col="blue")
title("Boxplots of logCPMs")

top_500_classic_F <- dge_output_classic_F[["table"]]
top500_lcpm_classic_F <- logcounts_F[na.omit(match(rownames(top_500_classic_F),rownames(logcounts_F))),]
top500_lcpm.df_classic_F <- as.data.frame(top500_lcpm_classic_F)
write.csv(top500_lcpm.df_classic_F, "~/Documents/RNASeqBH1-BH32/EdgeR/BH1-32T_top500_lcpm_classic_F.csv")

#PCA top 500 F classic----
pcadata_500_F<- pca(t(top500_lcpm_classic_F), ncomp=21)

#Calculate variance
eigs <- pcadata_500_F$sdev^2
PCA1_500_var <- eigs[1] / sum(eigs)
PCA1_500_var<-round(PCA1_500_var,digits = 2)*100
PCA2_500_var <- eigs[2] / sum(eigs)
PCA2_500_var<-round(PCA2_500_var,digits = 2)*100

pcadata_plot <- data.frame(pcadata_500_F$variates$X[,1:2])

ggplot(pcadata_plot,aes(PC1,PC2,colour=Y, shape=Y))+
  geom_point(size=7)+
  coord_fixed(ratio=1)+
  theme_classic(base_size=30)+
  labs(title= "e) PCA, top 500, classic, F", colour="Pubertal Status", shape="Pubertal Status")+
  xlab(label= paste("PC1 (expl. variance= ", PCA1_500_var, "%)",sep=""))+
  ylab(label= paste("PC2 (expl. variance= ", PCA2_500_var, "%)",sep=""))+
  theme(legend.position = "bottom")+
  scale_color_manual(values = c("#1F77B4", "#FF7F0E"))+
  stat_stars(aes(colour=Y))+
  stat_ellipse((aes(colour=Y)), level=0.95, linewidth = 1.5, type="t")
ggsave("SQePCAtop500genesFpubertyclassicNTA.png", width = 10, height = 10, dpi=320)

#Explained variance
plot(pcadata_500_F)
plotLoadings(pcadata_500_F)

#PLSDA top 500 F classic----
plsda_500_F<- plsda(t(top500_lcpm_classic_F), Y, ncomp = 21)

#Calculate variance
ev <- explained_variance(plsda_500_F$X, plsda_500_F$variates$X, ncomp =21)
Comp1v<-round(ev[1],digits = 2)*100
Comp2v<-round(ev[2],digits = 2)*100

plsdadata_plot <- data.frame(plsda_500_F$variates$X[,1:2])

ggplot(plsdadata_plot,aes(comp1,comp2,colour=Y, shape=Y))+
  geom_point(size=7)+
  coord_fixed(ratio = 1)+
  theme_classic(base_size=30)+
  labs(title= "g) PLSDA, top 500, classic, F", colour="Pubertal Status", shape="Pubertal Status")+
  xlab(label= paste("Component 1 (expl. variance= ", Comp1v, "%)",sep=""))+
  ylab(label= paste("Component 2 (expl. variance= ", Comp2v, "%)",sep=""))+
  theme(legend.position= "bottom")+
  scale_color_manual(values = c("#1F77B4", "#FF7F0E"))+
  stat_stars(aes(colour=Y))+
  stat_ellipse((aes(colour=Y)), level=0.95, linewidth = 1.5, type="t")
ggsave("SQgPLSDAtop500genesFpubertyclassicNTA.png", width = 10, height = 10, dpi=320)
