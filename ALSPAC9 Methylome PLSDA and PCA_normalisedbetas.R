setwd("~/Documents/ALSPAC Methylome")

library(mixOmics)
library(rgl)
library(tidyverse)
library(ggpubr)

#Load data----
load("~/Documents/ALSPAC Methylome/Data/samplesheet/samplesheet.Robj")
load("~/Documents/ALSPAC Methylome/Data/B2norm.Robj")

#Remove duplicate rows as per ALSPAC methylome user guide
Dups<-which(samplesheet$duplicate.rm=="Remove")
samplesheet<-samplesheet[-Dups,]

#Change row name to participant ID plus rank (combine columns 2 & 3, no spaces), as for Clin Data file
rownames(samplesheet)<-paste(samplesheet[,2],samplesheet[,3],sep="")

#B2 all probes----
ClinData_methylomecoB2 <- read.csv("~/Documents/ALSPAC_Histograms/ClinData_methylomecoB2.csv")
rownames(ClinData_methylomecoB2) <- ClinData_methylomecoB2$X
ClinData_methylomecoB2 <- ClinData_methylomecoB2[,-c(1:2)]

#Subset of samplesheet for participants with time to B2
samplesheetB2 <- samplesheet[match(rownames(ClinData_methylomecoB2),rownames(samplesheet)),]

#Add quartile column to & fill as NA
ClinData_methylomecoB2$DmB2Quartile <-NA 

#Calculate quartiles of time to B2 & put in new object
Quartiles_timetoB2 <- quantile(ClinData_methylomecoB2$DmB2)

#Populate time to B2 quartile column with 1st-4th quartile
ClinData_methylomecoB2$DmB2Quartile<- ifelse(ClinData_methylomecoB2$DmB2<Quartiles_timetoB2[2],"1st Quartile", ClinData_methylomecoB2$DmB2Quartile)
ClinData_methylomecoB2$DmB2Quartile<- ifelse(ClinData_methylomecoB2$DmB2>=Quartiles_timetoB2[2] & ClinData_methylomecoB2$DmB2<Quartiles_timetoB2[3], "2nd Quartile", ClinData_methylomecoB2$DmB2Quartile)
ClinData_methylomecoB2$DmB2Quartile<- ifelse(ClinData_methylomecoB2$DmB2>=Quartiles_timetoB2[3] & ClinData_methylomecoB2$DmB2<Quartiles_timetoB2[4], "3rd Quartile", ClinData_methylomecoB2$DmB2Quartile)
ClinData_methylomecoB2$DmB2Quartile<- ifelse(ClinData_methylomecoB2$DmB2>=Quartiles_timetoB2[4], "4th Quartile", ClinData_methylomecoB2$DmB2Quartile)

#Create a puberty groups file for RF
B2QMGp <- ClinData_methylomecoB2[,c(2290,2292)]
B2QMGp$ClQGp[B2QMGp$DmB2Quartile=="1st Quartile"] <- "1st Quartile"
B2QMGp$ClQGp[B2QMGp$DmB2Quartile=="2nd Quartile"] <- "2nd-4th Quartiles"
B2QMGp$ClQGp[B2QMGp$DmB2Quartile=="3rd Quartile"] <- "2nd-4th Quartiles"
B2QMGp$ClQGp[B2QMGp$DmB2Quartile=="4th Quartile"] <- "2nd-4th Quartiles"

B2QMGp$FuQGp[B2QMGp$DmB2Quartile=="1st Quartile"] <- "1st-3rd Quartiles"
B2QMGp$FuQGp[B2QMGp$DmB2Quartile=="2nd Quartile"] <- "1st-3rd Quartiles"
B2QMGp$FuQGp[B2QMGp$DmB2Quartile=="3rd Quartile"] <- "1st-3rd Quartiles"
B2QMGp$FuQGp[B2QMGp$DmB2Quartile=="4th Quartile"] <- "4th Quartile"
write.csv(B2QMGp,"B2QMGp.csv",row.names = TRUE)

#Convert quartile column in puberty groups object to a factor
Y <- as.factor(B2QMGp$DmB2Quartile)

#PLSDA, all probes, time to B2
plsda.dataY<- plsda(B2norm, Y, ncomp = 354)

#Calculate variance
ev <- explained_variance(plsda.dataY$X, plsda.dataY$variates$X, ncomp =354)
Comp1v<-round(ev[1],digits = 2)*100
Comp2v<-round(ev[2],digits = 2)*100

plsdadata_plot <- data.frame(plsda.dataY$variates$X[,1:2])

ggplot(plsdadata_plot,aes(comp1,comp2,colour=Y, shape=Y))+
  geom_point(size=5, alpha=0.8)+
  coord_fixed(ratio=1)+
  theme_classic(base_size=30)+
  labs(title= "b) PLSDA, all probes", colour="Time to B2", shape="Time to B2")+
  xlab(label= paste("Component 1 (expl. variance= ", Comp1v, "%)",sep=""))+
  ylab(label= paste("Component 2 (expl. variance= ", Comp2v, "%)",sep=""))+
  theme(legend.position= "bottom")+
  scale_color_manual(values = c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728"), labels=c("Q1", "Q2", "Q3", "Q4"))+
  scale_shape_manual(values=c(1, 2, 3, 4), labels=c("Q1", "Q2", "Q3", "Q4"))+
  stat_stars(aes(colour=Y), alpha=0.5)+
  stat_ellipse((aes(colour=Y)), level=0.95, linewidth = 1.5, type="t", alpha=0.5)
ggsave("SQbPLSDAallB2.png", width = 10, height = 10, dpi=320)

#PCA all probes, time to B2
Y.pca <- pca(B2norm, ncomp=354) 

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
  scale_color_manual(values = c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728"), labels=c("Q1", "Q2", "Q3", "Q4"))+
  scale_shape_manual(values=c(1, 2, 3, 4), labels=c("Q1", "Q2", "Q3", "Q4"))+
  stat_stars(aes(colour=Y), alpha=0.5)+
  stat_ellipse((aes(colour=Y)), level=0.95, linewidth = 1.5, type="t", alpha=0.5)
ggsave("SQaPCAallB2.png", width = 10, height = 10, dpi=320)

#B2, top 500 probes----

#Closest to B2
Top500B2_Cl <- read.csv("~/Documents/ALSPAC Methylome/Data/DMPcB2top500.csv")
Top500B2_Cl <- Top500B2_Cl[,1]

#Subset betasB2 for top 500_Cl
betasB2top500_Cl <- B2norm[,match(Top500B2_Cl,colnames(B2norm))]
write.csv(betasB2top500_Cl,"Top500B2_Cl_norm.csv")

#PLSDA, top 500 probes closest to B2
plsda.dataY<- plsda(betasB2top500_Cl, Y, ncomp = 354)

#Calculate variance
ev <- explained_variance(plsda.dataY$X, plsda.dataY$variates$X, ncomp =354)
Comp1v<-round(ev[1],digits = 2)*100
Comp2v<-round(ev[2],digits = 2)*100

plsdadata_plot <- data.frame(plsda.dataY$variates$X[,1:2])

ggplot(plsdadata_plot,aes(comp1,comp2,colour=Y, shape=Y))+
  geom_point(size=5, alpha=0.8)+
  coord_fixed(ratio=1)+
  theme_classic(base_size=30)+
  labs(title= "e) PLSDA, top 500, closest", colour="Time to B2", shape="Time to B2")+
  xlab(label= paste("Component 1 (expl. variance= ", Comp1v, "%)",sep=""))+
  ylab(label= paste("Component 2 (expl. variance= ", Comp2v, "%)",sep=""))+
  theme(legend.position= "bottom")+
  scale_color_manual(values = c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728"), labels=c("Q1", "Q2", "Q3", "Q4"))+
  scale_shape_manual(values=c(1, 2, 3, 4), labels=c("Q1", "Q2", "Q3", "Q4"))+
  stat_stars(aes(colour=Y), alpha=0.5)+
  stat_ellipse((aes(colour=Y)), level=0.95, linewidth = 1.5, type="t", alpha=0.5)
ggsave("SQePLSDA500clB2.png", width = 10, height = 10, dpi=320)

#PCA, top 500 probes closest to B2
Y.pca <- pca(betasB2top500_Cl, ncomp=354) 

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
  labs(title= "c) PCA, top 500, closest", colour="Time to B2", shape="Time to B2")+
  xlab(label= paste("PC1 (expl. variance= ", PCA1_var, "%)",sep=""))+
  ylab(label= paste("PC2 (expl. variance= ", PCA2_var, "%)",sep=""))+
  theme(legend.position = "bottom")+
  scale_color_manual(values = c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728"), labels=c("Q1", "Q2", "Q3", "Q4"))+
  scale_shape_manual(values=c(1, 2, 3, 4), labels=c("Q1", "Q2", "Q3", "Q4"))+
  stat_stars(aes(colour=Y), alpha=0.5)+
  stat_ellipse((aes(colour=Y)), level=0.95, linewidth = 1.5, type="t", alpha=0.5)
ggsave("SQcPCA500clB2.png", width = 10, height = 10, dpi=320)

#Furthest from B2
Top500B2_Fu <- read.csv("~/Documents/ALSPAC Methylome/Data/DMPfB2top500.csv")
Top500B2_Fu <- Top500B2_Fu[,1]

#Subset betasB2 for top 500_Fu
betasB2top500_Fu <- B2norm[,match(Top500B2_Fu,colnames(B2norm))]
write.csv(betasB2top500_Fu,"Top500B2_Fu_norm.csv")

#PLSDA, top 500 probes furthest from B2
plsda.dataY<- plsda(betasB2top500_Fu, Y, ncomp = 354)

#Calculate variance
ev <- explained_variance(plsda.dataY$X, plsda.dataY$variates$X, ncomp =354)
Comp1v<-round(ev[1],digits = 2)*100
Comp2v<-round(ev[2],digits = 2)*100

plsdadata_plot <- data.frame(plsda.dataY$variates$X[,1:2])

ggplot(plsdadata_plot,aes(comp1,comp2,colour=Y, shape=Y))+
  geom_point(size=5, alpha=0.8)+
  coord_fixed(ratio=1)+
  theme_classic(base_size=30)+
  labs(title= "f) PLSDA, top 500, furthest", colour="Time to B2", shape="Time to B2")+
  xlab(label= paste("Component 1 (expl. variance= ", Comp1v, "%)",sep=""))+
  ylab(label= paste("Component 2 (expl. variance= ", Comp2v, "%)",sep=""))+
  theme(legend.position= "bottom")+
  scale_color_manual(values = c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728"), labels=c("Q1", "Q2", "Q3", "Q4"))+
  scale_shape_manual(values=c(1, 2, 3, 4), labels=c("Q1", "Q2", "Q3", "Q4"))+
  stat_stars(aes(colour=Y), alpha=0.5)+
  stat_ellipse((aes(colour=Y)), level=0.95, linewidth = 1.5, type="t", alpha=0.5)
ggsave("SQfPLSDA500fuB2.png", width = 10, height = 10, dpi=320)

#PCA, top 500 probes furthest from B2
Y.pca <- pca(betasB2top500_Fu, ncomp=354) 

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
  labs(title= "d) PCA, top 500, furthest", colour="Time to B2", shape="Time to B2")+
  xlab(label= paste("PC1 (expl. variance= ", PCA1_var, "%)",sep=""))+
  ylab(label= paste("PC2 (expl. variance= ", PCA2_var, "%)",sep=""))+
  theme(legend.position = "bottom")+
  scale_color_manual(values = c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728"), labels=c("Q1", "Q2", "Q3", "Q4"))+
  scale_shape_manual(values=c(1, 2, 3, 4), labels=c("Q1", "Q2", "Q3", "Q4"))+
  stat_stars(aes(colour=Y), alpha=0.5)+
  stat_ellipse((aes(colour=Y)), level=0.95, linewidth = 1.5, type="t", alpha=0.5)
ggsave("SQdPCA500fuB2.png", width = 10, height = 10, dpi=320)

rm(B2QMGp)
rm(B2norm)
rm(betasB2top500_Cl)
rm(betasB2top500_Fu)
rm(ClinData_methylomecoB2)
rm(samplesheetB2)

#P2M all probes----
load("~/Documents/ALSPAC Methylome/Data/P2Mnormc.Robj")
ClinData_methylomecoP2M <- read.csv("~/Documents/ALSPAC_Histograms/ClinData_methylomecoP2M.csv")
rownames(ClinData_methylomecoP2M) <- ClinData_methylomecoP2M$X
ClinData_methylomecoP2M <- ClinData_methylomecoP2M[,-c(1:2)]

#Subset of samplesheet for participants with time to P2M
samplesheetP2M <- samplesheet[match(rownames(ClinData_methylomecoP2M),rownames(samplesheet)),]

#Add quartile column to & fill as NA
ClinData_methylomecoP2M$DmP2MQuartile <-NA 

#Calculate quartiles of time to P2M & put in new object
Quartiles_timetoP2M <- quantile(ClinData_methylomecoP2M$DmP2M)

#Populate time to P2M quartile column with 1st-4th quartile
ClinData_methylomecoP2M$DmP2MQuartile<- ifelse(ClinData_methylomecoP2M$DmP2M<Quartiles_timetoP2M[2],"1st Quartile", ClinData_methylomecoP2M$DmP2MQuartile)
ClinData_methylomecoP2M$DmP2MQuartile<- ifelse(ClinData_methylomecoP2M$DmP2M>=Quartiles_timetoP2M[2] & ClinData_methylomecoP2M$DmP2M<Quartiles_timetoP2M[3], "2nd Quartile", ClinData_methylomecoP2M$DmP2MQuartile)
ClinData_methylomecoP2M$DmP2MQuartile<- ifelse(ClinData_methylomecoP2M$DmP2M>=Quartiles_timetoP2M[3] & ClinData_methylomecoP2M$DmP2M<Quartiles_timetoP2M[4], "3rd Quartile", ClinData_methylomecoP2M$DmP2MQuartile)
ClinData_methylomecoP2M$DmP2MQuartile<- ifelse(ClinData_methylomecoP2M$DmP2M>=Quartiles_timetoP2M[4], "4th Quartile", ClinData_methylomecoP2M$DmP2MQuartile)

#Create a puberty groups file for RF
P2MQMGp <- ClinData_methylomecoP2M[,c(2288,2292)]
P2MQMGp$ClQGp[P2MQMGp$DmP2MQuartile=="1st Quartile"] <- "1st Quartile"
P2MQMGp$ClQGp[P2MQMGp$DmP2MQuartile=="2nd Quartile"] <- "2nd-4th Quartiles"
P2MQMGp$ClQGp[P2MQMGp$DmP2MQuartile=="3rd Quartile"] <- "2nd-4th Quartiles"
P2MQMGp$ClQGp[P2MQMGp$DmP2MQuartile=="4th Quartile"] <- "2nd-4th Quartiles"

P2MQMGp$FuQGp[P2MQMGp$DmP2MQuartile=="1st Quartile"] <- "1st-3rd Quartiles"
P2MQMGp$FuQGp[P2MQMGp$DmP2MQuartile=="2nd Quartile"] <- "1st-3rd Quartiles"
P2MQMGp$FuQGp[P2MQMGp$DmP2MQuartile=="3rd Quartile"] <- "1st-3rd Quartiles"
P2MQMGp$FuQGp[P2MQMGp$DmP2MQuartile=="4th Quartile"] <- "4th Quartile"
write.csv(P2MQMGp,"P2MQMGp.csv",row.names = TRUE)

#Convert quartile column in puberty groups object to a factor
Y <- as.factor(P2MQMGp$DmP2MQuartile)

#PLSDA, all probes, time to P2M
plsda.dataY<- plsda(P2Mnormc, Y, ncomp=242)

#Calculate variance
ev <- explained_variance(plsda.dataY$X, plsda.dataY$variates$X, ncomp=242)
Comp1v<-round(ev[1],digits = 2)*100
Comp2v<-round(ev[2],digits = 2)*100

plsdadata_plot <- data.frame(plsda.dataY$variates$X[,1:2])

ggplot(plsdadata_plot,aes(comp1,comp2,colour=Y, shape=Y))+
  geom_point(size=5, alpha=0.8)+
  coord_fixed(ratio=1)+
  theme_classic(base_size=30)+
  labs(title= "b) PLSDA, all probes", colour="Time to P2M", shape="Time to P2M")+
  xlab(label= paste("Component 1 (expl. variance= ", Comp1v, "%)",sep=""))+
  ylab(label= paste("Component 2 (expl. variance= ", Comp2v, "%)",sep=""))+
  theme(legend.position= "bottom")+
  scale_color_manual(values = c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728"), labels=c("Q1", "Q2", "Q3", "Q4"))+
  scale_shape_manual(values=c(1, 2, 3, 4), labels=c("Q1", "Q2", "Q3", "Q4"))+
  stat_stars(aes(colour=Y), alpha=0.5)+
  stat_ellipse((aes(colour=Y)), level=0.95, linewidth = 1.5, type="t", alpha=0.5)
ggsave("SQbPLSDAallP2M.png", width = 10, height = 10, dpi=320)

#PCA all probes, time to P2M
Y.pca <- pca(P2Mnormc, ncomp=242) 

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
  labs(title= "a) PCA, all probes", colour="Time to P2M", shape="Time to P2M")+
  xlab(label= paste("PC1 (expl. variance= ", PCA1_var, "%)",sep=""))+
  ylab(label= paste("PC2 (expl. variance= ", PCA2_var, "%)",sep=""))+
  theme(legend.position = "bottom")+
  scale_color_manual(values = c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728"), labels=c("Q1", "Q2", "Q3", "Q4"))+
  scale_shape_manual(values=c(1, 2, 3, 4), labels=c("Q1", "Q2", "Q3", "Q4"))+
  stat_stars(aes(colour=Y), alpha=0.5)+
  stat_ellipse((aes(colour=Y)), level=0.95, linewidth = 1.5, type="t", alpha=0.5)
ggsave("SQaPCAallP2M.png", width = 10, height = 10, dpi=320)

#P2M, top 500 probes----

#Closest to P2M
Top500P2M_Cl <- read.csv("~/Documents/ALSPAC Methylome/Data/DMPcP2Mtop500.csv")
Top500P2M_Cl <- Top500P2M_Cl[,1]

#Subset betasP2M for top 500_Cl
betasP2Mtop500_Cl <- P2Mnormc[,match(Top500P2M_Cl,colnames(P2Mnormc))]
write.csv(betasP2Mtop500_Cl,"Top500P2M_Cl_norm.csv")

#PLSDA, top 500 probes closest to P2M
plsda.dataY<- plsda(betasP2Mtop500_Cl, Y, ncomp=242)

#Calculate variance
ev <- explained_variance(plsda.dataY$X, plsda.dataY$variates$X, ncomp=242)
Comp1v<-round(ev[1],digits = 2)*100
Comp2v<-round(ev[2],digits = 2)*100

plsdadata_plot <- data.frame(plsda.dataY$variates$X[,1:2])

ggplot(plsdadata_plot,aes(comp1,comp2,colour=Y, shape=Y))+
  geom_point(size=5, alpha=0.8)+
  coord_fixed(ratio=1)+
  theme_classic(base_size=30)+
  labs(title= "e) PLSDA, top 500, closest", colour="Time to P2M", shape="Time to P2M")+
  xlab(label= paste("Component 1 (expl. variance= ", Comp1v, "%)",sep=""))+
  ylab(label= paste("Component 2 (expl. variance= ", Comp2v, "%)",sep=""))+
  theme(legend.position= "bottom")+
  scale_color_manual(values = c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728"), labels=c("Q1", "Q2", "Q3", "Q4"))+
  scale_shape_manual(values=c(1, 2, 3, 4), labels=c("Q1", "Q2", "Q3", "Q4"))+
  stat_stars(aes(colour=Y), alpha=0.5)+
  stat_ellipse((aes(colour=Y)), level=0.95, linewidth = 1.5, type="t", alpha=0.5)
ggsave("SQePLSDA500clP2M.png", width = 10, height = 10, dpi=320)

#PCA, top 500 probes closest to P2M
Y.pca <- pca(betasP2Mtop500_Cl, ncomp=242) 

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
  labs(title= "c) PCA, top 500, closest", colour="Time to P2M", shape="Time to P2M")+
  xlab(label= paste("PC1 (expl. variance= ", PCA1_var, "%)",sep=""))+
  ylab(label= paste("PC2 (expl. variance= ", PCA2_var, "%)",sep=""))+
  theme(legend.position = "bottom")+
  scale_color_manual(values = c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728"), labels=c("Q1", "Q2", "Q3", "Q4"))+
  scale_shape_manual(values=c(1, 2, 3, 4), labels=c("Q1", "Q2", "Q3", "Q4"))+
  stat_stars(aes(colour=Y), alpha=0.5)+
  stat_ellipse((aes(colour=Y)), level=0.95, linewidth = 1.5, type="t", alpha=0.5)
ggsave("SQcPCA500clP2M.png", width = 10, height = 10, dpi=320)

#Furthest from P2M
Top500P2M_Fu <- read.csv("~/Documents/ALSPAC Methylome/Data/DMPfP2Mtop500.csv")
Top500P2M_Fu <- Top500P2M_Fu[,1]

#Subset betasP2M for top 500_Fu
betasP2Mtop500_Fu <- P2Mnormc[,match(Top500P2M_Fu,colnames(P2Mnormc))]
write.csv(betasP2Mtop500_Fu,"Top500P2M_Fu_norm.csv")

#PLSDA, top 500 probes furthest from P2M
plsda.dataY<- plsda(betasP2Mtop500_Fu, Y, ncomp=242)

#Calculate variance
ev <- explained_variance(plsda.dataY$X, plsda.dataY$variates$X, ncomp=242)
Comp1v<-round(ev[1],digits = 2)*100
Comp2v<-round(ev[2],digits = 2)*100

plsdadata_plot <- data.frame(plsda.dataY$variates$X[,1:2])

ggplot(plsdadata_plot,aes(comp1,comp2,colour=Y, shape=Y))+
  geom_point(size=5, alpha=0.8)+
  coord_fixed(ratio=1)+
  theme_classic(base_size=30)+
  labs(title= "f) PLSDA, top 500, furthest", colour="Time to P2M", shape="Time to P2M")+
  xlab(label= paste("Component 1 (expl. variance= ", Comp1v, "%)",sep=""))+
  ylab(label= paste("Component 2 (expl. variance= ", Comp2v, "%)",sep=""))+
  theme(legend.position= "bottom")+
  scale_color_manual(values = c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728"), labels=c("Q1", "Q2", "Q3", "Q4"))+
  scale_shape_manual(values=c(1, 2, 3, 4), labels=c("Q1", "Q2", "Q3", "Q4"))+
  stat_stars(aes(colour=Y), alpha=0.5)+
  stat_ellipse((aes(colour=Y)), level=0.95, linewidth = 1.5, type="t", alpha=0.5)
ggsave("SQfPLSDA500fuP2M.png", width = 10, height = 10, dpi=320)

#PCA, top 500 probes furthest from P2M
Y.pca <- pca(betasP2Mtop500_Fu, ncomp=242)

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
  labs(title= "d) PCA, top 500, furthest", colour="Time to P2M", shape="Time to P2M")+
  xlab(label= paste("PC1 (expl. variance= ", PCA1_var, "%)",sep=""))+
  ylab(label= paste("PC2 (expl. variance= ", PCA2_var, "%)",sep=""))+
  theme(legend.position = "bottom")+
  scale_color_manual(values = c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728"), labels=c("Q1", "Q2", "Q3", "Q4"))+
  scale_shape_manual(values=c(1, 2, 3, 4), labels=c("Q1", "Q2", "Q3", "Q4"))+
  stat_stars(aes(colour=Y), alpha=0.5)+
  stat_ellipse((aes(colour=Y)), level=0.95, linewidth = 1.5, type="t", alpha=0.5)
ggsave("SQdPCA500fuP2M.png", width = 10, height = 10, dpi=320)

rm(P2MQMGp)
rm(P2Mnormc)
rm(betasP2Mtop500_Cl)
rm(betasP2Mtop500_Fu)
rm(ClinData_methylomecoP2M)
rm(samplesheetP2M)

#P2F all probes----
load("~/Documents/ALSPAC Methylome/Data/P2Fnormc.Robj")
ClinData_methylomecoP2F <- read.csv("~/Documents/ALSPAC_Histograms/ClinData_methylomecoP2F.csv")
rownames(ClinData_methylomecoP2F) <- ClinData_methylomecoP2F$X
ClinData_methylomecoP2F <- ClinData_methylomecoP2F[,-c(1:2)]

#Subset of samplesheet for participants with time to P2F
samplesheetP2F <- samplesheet[match(rownames(ClinData_methylomecoP2F),rownames(samplesheet)),]

#Add quartile column to & fill as NA
ClinData_methylomecoP2F$DmP2FQuartile <-NA 

#Calculate quartiles of time to P2F & put in new object
Quartiles_timetoP2F <- quantile(ClinData_methylomecoP2F$DmP2F)

#Populate time to P2F quartile column with 1st-4th quartile
ClinData_methylomecoP2F$DmP2FQuartile<- ifelse(ClinData_methylomecoP2F$DmP2F<Quartiles_timetoP2F[2],"1st Quartile", ClinData_methylomecoP2F$DmP2FQuartile)
ClinData_methylomecoP2F$DmP2FQuartile<- ifelse(ClinData_methylomecoP2F$DmP2F>=Quartiles_timetoP2F[2] & ClinData_methylomecoP2F$DmP2F<Quartiles_timetoP2F[3], "2nd Quartile", ClinData_methylomecoP2F$DmP2FQuartile)
ClinData_methylomecoP2F$DmP2FQuartile<- ifelse(ClinData_methylomecoP2F$DmP2F>=Quartiles_timetoP2F[3] & ClinData_methylomecoP2F$DmP2F<Quartiles_timetoP2F[4], "3rd Quartile", ClinData_methylomecoP2F$DmP2FQuartile)
ClinData_methylomecoP2F$DmP2FQuartile<- ifelse(ClinData_methylomecoP2F$DmP2F>=Quartiles_timetoP2F[4], "4th Quartile", ClinData_methylomecoP2F$DmP2FQuartile)

#Create a puberty groups file for RF
P2FQMGp <- ClinData_methylomecoP2F[,c(2289,2292)]
P2FQMGp$ClQGp[P2FQMGp$DmP2FQuartile=="1st Quartile"] <- "1st Quartile"
P2FQMGp$ClQGp[P2FQMGp$DmP2FQuartile=="2nd Quartile"] <- "2nd-4th Quartiles"
P2FQMGp$ClQGp[P2FQMGp$DmP2FQuartile=="3rd Quartile"] <- "2nd-4th Quartiles"
P2FQMGp$ClQGp[P2FQMGp$DmP2FQuartile=="4th Quartile"] <- "2nd-4th Quartiles"

P2FQMGp$FuQGp[P2FQMGp$DmP2FQuartile=="1st Quartile"] <- "1st-3rd Quartiles"
P2FQMGp$FuQGp[P2FQMGp$DmP2FQuartile=="2nd Quartile"] <- "1st-3rd Quartiles"
P2FQMGp$FuQGp[P2FQMGp$DmP2FQuartile=="3rd Quartile"] <- "1st-3rd Quartiles"
P2FQMGp$FuQGp[P2FQMGp$DmP2FQuartile=="4th Quartile"] <- "4th Quartile"
write.csv(P2FQMGp,"P2FQMGp.csv",row.names = TRUE)

#Convert quartile column in puberty groups object to a factor
Y <- as.factor(P2FQMGp$DmP2FQuartile)

#PLSDA, all probes, time to P2F
plsda.dataY<- plsda(P2Fnormc, Y, ncomp=273)

#Calculate variance
ev <- explained_variance(plsda.dataY$X, plsda.dataY$variates$X, ncomp=273)
Comp1v<-round(ev[1],digits = 2)*100
Comp2v<-round(ev[2],digits = 2)*100

plsdadata_plot <- data.frame(plsda.dataY$variates$X[,1:2])

ggplot(plsdadata_plot,aes(comp1,comp2,colour=Y, shape=Y))+
  geom_point(size=5, alpha=0.8)+
  coord_fixed(ratio=1)+
  theme_classic(base_size=30)+
  labs(title= "b) PLSDA, all probes", colour="Time to P2F", shape="Time to P2F")+
  xlab(label= paste("Component 1 (expl. variance= ", Comp1v, "%)",sep=""))+
  ylab(label= paste("Component 2 (expl. variance= ", Comp2v, "%)",sep=""))+
  theme(legend.position= "bottom")+
  scale_color_manual(values = c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728"), labels=c("Q1", "Q2", "Q3", "Q4"))+
  scale_shape_manual(values=c(1, 2, 3, 4), labels=c("Q1", "Q2", "Q3", "Q4"))+
  stat_stars(aes(colour=Y), alpha=0.5)+
  stat_ellipse((aes(colour=Y)), level=0.95, linewidth = 1.5, type="t", alpha=0.5)
ggsave("SQbPLSDAallP2F.png", width = 10, height = 10, dpi=320)

#PCA all probes, time to P2F
Y.pca <- pca(P2Fnormc, ncomp=273) 

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
  labs(title= "a) PCA, all probes", colour="Time to P2F", shape="Time to P2F")+
  xlab(label= paste("PC1 (expl. variance= ", PCA1_var, "%)",sep=""))+
  ylab(label= paste("PC2 (expl. variance= ", PCA2_var, "%)",sep=""))+
  theme(legend.position = "bottom")+
  scale_color_manual(values = c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728"), labels=c("Q1", "Q2", "Q3", "Q4"))+
  scale_shape_manual(values=c(1, 2, 3, 4), labels=c("Q1", "Q2", "Q3", "Q4"))+
  stat_stars(aes(colour=Y), alpha=0.5)+
  stat_ellipse((aes(colour=Y)), level=0.95, linewidth = 1.5, type="t", alpha=0.5)
ggsave("SQaPCAallP2F.png", width = 10, height = 10, dpi=320)

#P2F, top 500 probes----

#Closest to P2F
Top500P2F_Cl <- read.csv("~/Documents/ALSPAC Methylome/Data/DMPcP2Ftop500.csv")
Top500P2F_Cl <- Top500P2F_Cl[,1]

#Subset betasP2F for top 500_Cl
betasP2Ftop500_Cl <- P2Fnormc[,match(Top500P2F_Cl,colnames(P2Fnormc))]
write.csv(betasP2Ftop500_Cl,"Top500P2F_Cl_norm.csv")

#PLSDA, top 500 probes closest to P2F
plsda.dataY<- plsda(betasP2Ftop500_Cl, Y, ncomp=273)

#Calculate variance
ev <- explained_variance(plsda.dataY$X, plsda.dataY$variates$X, ncomp=273)
Comp1v<-round(ev[1],digits = 2)*100
Comp2v<-round(ev[2],digits = 2)*100

plsdadata_plot <- data.frame(plsda.dataY$variates$X[,1:2])

ggplot(plsdadata_plot,aes(comp1,comp2,colour=Y, shape=Y))+
  geom_point(size=5, alpha=0.8)+
  coord_fixed(ratio=1)+
  theme_classic(base_size=30)+
  labs(title= "e) PLSDA, top 500, closest", colour="Time to P2F", shape="Time to P2F")+
  xlab(label= paste("Component 1 (expl. variance= ", Comp1v, "%)",sep=""))+
  ylab(label= paste("Component 2 (expl. variance= ", Comp2v, "%)",sep=""))+
  theme(legend.position= "bottom")+
  scale_color_manual(values = c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728"), labels=c("Q1", "Q2", "Q3", "Q4"))+
  scale_shape_manual(values=c(1, 2, 3, 4), labels=c("Q1", "Q2", "Q3", "Q4"))+
  stat_stars(aes(colour=Y), alpha=0.5)+
  stat_ellipse((aes(colour=Y)), level=0.95, linewidth = 1.5, type="t", alpha=0.5)
ggsave("SQePLSDA500clP2F.png", width = 10, height = 10, dpi=320)

#PCA, top 500 probes closest to P2F
Y.pca <- pca(betasP2Ftop500_Cl, ncomp=273) 

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
  labs(title= "c) PCA, top 500, closest", colour="Time to P2F", shape="Time to P2F")+
  xlab(label= paste("PC1 (expl. variance= ", PCA1_var, "%)",sep=""))+
  ylab(label= paste("PC2 (expl. variance= ", PCA2_var, "%)",sep=""))+
  theme(legend.position = "bottom")+
  scale_color_manual(values = c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728"), labels=c("Q1", "Q2", "Q3", "Q4"))+
  scale_shape_manual(values=c(1, 2, 3, 4), labels=c("Q1", "Q2", "Q3", "Q4"))+
  stat_stars(aes(colour=Y), alpha=0.5)+
  stat_ellipse((aes(colour=Y)), level=0.95, linewidth = 1.5, type="t", alpha=0.5)
ggsave("SQcPCA500clP2F.png", width = 10, height = 10, dpi=320)

#Furthest from P2F
Top500P2F_Fu <- read.csv("~/Documents/ALSPAC Methylome/Data/DMPfP2Ftop500.csv")
Top500P2F_Fu <- Top500P2F_Fu[,1]

#Subset betasP2F for top 500_Fu
betasP2Ftop500_Fu <- P2Fnormc[,match(Top500P2F_Fu,colnames(P2Fnormc))]
write.csv(betasP2Ftop500_Fu,"Top500P2F_Fu_norm.csv")

#PLSDA, top 500 probes furthest from P2F
plsda.dataY<- plsda(betasP2Ftop500_Fu, Y, ncomp=273)

#Calculate variance
ev <- explained_variance(plsda.dataY$X, plsda.dataY$variates$X, ncomp=273)
Comp1v<-round(ev[1],digits = 2)*100
Comp2v<-round(ev[2],digits = 2)*100

plsdadata_plot <- data.frame(plsda.dataY$variates$X[,1:2])

ggplot(plsdadata_plot,aes(comp1,comp2,colour=Y, shape=Y))+
  geom_point(size=5, alpha=0.8)+
  coord_fixed(ratio=1)+
  theme_classic(base_size=30)+
  labs(title= "f) PLSDA, top 500, furthest", colour="Time to P2F", shape="Time to P2F")+
  xlab(label= paste("Component 1 (expl. variance= ", Comp1v, "%)",sep=""))+
  ylab(label= paste("Component 2 (expl. variance= ", Comp2v, "%)",sep=""))+
  theme(legend.position= "bottom")+
  scale_color_manual(values = c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728"), labels=c("Q1", "Q2", "Q3", "Q4"))+
  scale_shape_manual(values=c(1, 2, 3, 4), labels=c("Q1", "Q2", "Q3", "Q4"))+
  stat_stars(aes(colour=Y), alpha=0.5)+
  stat_ellipse((aes(colour=Y)), level=0.95, linewidth = 1.5, type="t", alpha=0.5)
ggsave("SQfPLSDA500fuP2F.png", width = 10, height = 10, dpi=320)

#PCA, top 500 probes furthest from P2F
Y.pca <- pca(betasP2Ftop500_Fu, ncomp=273) 

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
  labs(title= "d) PCA, top 500, furthest", colour="Time to P2F", shape="Time to P2F")+
  xlab(label= paste("PC1 (expl. variance= ", PCA1_var, "%)",sep=""))+
  ylab(label= paste("PC2 (expl. variance= ", PCA2_var, "%)",sep=""))+
  theme(legend.position = "bottom")+
  scale_color_manual(values = c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728"), labels=c("Q1", "Q2", "Q3", "Q4"))+
  scale_shape_manual(values=c(1, 2, 3, 4), labels=c("Q1", "Q2", "Q3", "Q4"))+
  stat_stars(aes(colour=Y), alpha=0.5)+
  stat_ellipse((aes(colour=Y)), level=0.95, linewidth = 1.5, type="t", alpha=0.5)
ggsave("SQdPCA500fuP2F.png", width = 10, height = 10, dpi=320)

rm(P2FQMGp)
rm(P2Fnormc)
rm(betasP2Ftop500_Cl)
rm(betasP2Ftop500_Fu)
rm(ClinData_methylomecoP2F)
rm(samplesheetP2F)

#Menarche all probes----
load("~/Documents/ALSPAC Methylome/Data/Mennormc.Robj")
ClinData_methylomecoMen <- read.csv("~/Documents/Files from PEC drive/ALSPAC_Histograms/ClinData_methylomecoMen.csv")
rownames(ClinData_methylomecoMen) <- ClinData_methylomecoMen$X
ClinData_methylomecoMen <- ClinData_methylomecoMen[,-c(1:2)]

#Subset of samplesheet for participants with time to menarche
SamplesheetMen <- samplesheet[match(rownames(ClinData_methylomecoMen),rownames(samplesheet)),]

#Add quartile column to & fill as NA
ClinData_methylomecoMen$DmMenQuartile <-NA 

#Calculate quartiles of time to P2F & put in new object
Quartiles_timetoMen <- quantile(ClinData_methylomecoMen$DmMen)

#Populate time to menarche quartile column with 1st-4th quartile
ClinData_methylomecoMen$DmMenQuartile<- ifelse(ClinData_methylomecoMen$DmMen<Quartiles_timetoMen[2],"1st Quartile", ClinData_methylomecoMen$DmMenQuartile)
ClinData_methylomecoMen$DmMenQuartile<- ifelse(ClinData_methylomecoMen$DmMen>=Quartiles_timetoMen[2] & ClinData_methylomecoMen$DmMen<Quartiles_timetoMen[3], "2nd Quartile", ClinData_methylomecoMen$DmMenQuartile)
ClinData_methylomecoMen$DmMenQuartile<- ifelse(ClinData_methylomecoMen$DmMen>=Quartiles_timetoMen[3] & ClinData_methylomecoMen$DmMen<Quartiles_timetoMen[4], "3rd Quartile", ClinData_methylomecoMen$DmMenQuartile)
ClinData_methylomecoMen$DmMenQuartile<- ifelse(ClinData_methylomecoMen$DmMen>=Quartiles_timetoMen[4], "4th Quartile", ClinData_methylomecoMen$DmMenQuartile)

#Create a puberty groups file for RF
MenQMGp <- ClinData_methylomecoMen[,c(2291,2292)]
MenQMGp$ClQGp[MenQMGp$DmMenQuartile=="1st Quartile"] <- "1st Quartile"
MenQMGp$ClQGp[MenQMGp$DmMenQuartile=="2nd Quartile"] <- "2nd-4th Quartiles"
MenQMGp$ClQGp[MenQMGp$DmMenQuartile=="3rd Quartile"] <- "2nd-4th Quartiles"
MenQMGp$ClQGp[MenQMGp$DmMenQuartile=="4th Quartile"] <- "2nd-4th Quartiles"

MenQMGp$FuQGp[MenQMGp$DmMenQuartile=="1st Quartile"] <- "1st-3rd Quartiles"
MenQMGp$FuQGp[MenQMGp$DmMenQuartile=="2nd Quartile"] <- "1st-3rd Quartiles"
MenQMGp$FuQGp[MenQMGp$DmMenQuartile=="3rd Quartile"] <- "1st-3rd Quartiles"
MenQMGp$FuQGp[MenQMGp$DmMenQuartile=="4th Quartile"] <- "4th Quartile"
write.csv(MenQMGp,"MenQMGp.csv",row.names = TRUE)

#Convert quartile column in puberty groups object to a factor
Y <- as.factor(MenQMGp$DmMenQuartile)

#PLSDA, all probes, time to menarche
plsda.dataY<- plsda(Mennormc, Y, ncomp=463)

#Calculate variance
ev <- explained_variance(plsda.dataY$X, plsda.dataY$variates$X, ncomp=463)
Comp1v<-round(ev[1],digits = 2)*100
Comp2v<-round(ev[2],digits = 2)*100

plsdadata_plot <- data.frame(plsda.dataY$variates$X[,1:2])

ggplot(plsdadata_plot,aes(comp1,comp2,colour=Y, shape=Y))+
  geom_point(size=5, alpha=0.8)+
  coord_fixed(ratio=1)+
  theme_classic(base_size=30)+
  labs(title= "b) PLSDA, all probes", colour="Time to Menarche", shape="Time to Menarche")+
  xlab(label= paste("Component 1 (expl. variance= ", Comp1v, "%)",sep=""))+
  ylab(label= paste("Component 2 (expl. variance= ", Comp2v, "%)",sep=""))+
  theme(legend.position= "bottom")+
  scale_color_manual(values = c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728"), labels=c("Q1", "Q2", "Q3", "Q4"))+
  scale_shape_manual(values=c(1, 2, 3, 4), labels=c("Q1", "Q2", "Q3", "Q4"))+
  stat_stars(aes(colour=Y), alpha=0.5)+
  stat_ellipse((aes(colour=Y)), level=0.95, linewidth = 1.5, type="t", alpha=0.5)
ggsave("SQbPLSDAallMen.png", width = 10, height = 10, dpi=320)

#PCA all probes, time to menarche
Y.pca <- pca(Mennormc, ncomp=463) 

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
  labs(title= "a) PCA, all probes", colour="Time to Menarche", shape="Time to Menarche")+
  xlab(label= paste("PC1 (expl. variance= ", PCA1_var, "%)",sep=""))+
  ylab(label= paste("PC2 (expl. variance= ", PCA2_var, "%)",sep=""))+
  theme(legend.position = "bottom")+
  scale_color_manual(values = c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728"), labels=c("Q1", "Q2", "Q3", "Q4"))+
  scale_shape_manual(values=c(1, 2, 3, 4), labels=c("Q1", "Q2", "Q3", "Q4"))+
  stat_stars(aes(colour=Y), alpha=0.5)+
  stat_ellipse((aes(colour=Y)), level=0.95, linewidth = 1.5, type="t", alpha=0.5)
ggsave("SQaPCAallMen.png", width = 10, height = 10, dpi=320)

#Menarche, top 500 probes----

#Closest to menarche
Top500Men_Cl <- read.csv("~/Documents/ALSPAC Methylome/Data/DMPcMentop500.csv")
Top500Men_Cl <- Top500Men_Cl[,1]

#Subset betasMen for top 500_Cl
betasMentop500_Cl <- Mennormc[,match(Top500Men_Cl,colnames(Mennormc))]
write.csv(betasMentop500_Cl,"Top500Men_Cl_norm.csv")

#PLSDA, top 500 probes closest to menarche
plsda.dataY<- plsda(betasMentop500_Cl, Y, ncomp=463)

#Calculate variance
ev <- explained_variance(plsda.dataY$X, plsda.dataY$variates$X, ncomp=463)
Comp1v<-round(ev[1],digits = 2)*100
Comp2v<-round(ev[2],digits = 2)*100

plsdadata_plot <- data.frame(plsda.dataY$variates$X[,1:2])

ggplot(plsdadata_plot,aes(comp1,comp2,colour=Y, shape=Y))+
  geom_point(size=5, alpha=0.8)+
  coord_fixed(ratio=1)+
  theme_classic(base_size=30)+
  labs(title= "e) PLSDA, top 500, closest", colour="Time to Menarche", shape="Time to Menarche")+
  xlab(label= paste("Component 1 (expl. variance= ", Comp1v, "%)",sep=""))+
  ylab(label= paste("Component 2 (expl. variance= ", Comp2v, "%)",sep=""))+
  theme(legend.position= "bottom")+
  scale_color_manual(values = c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728"), labels=c("Q1", "Q2", "Q3", "Q4"))+
  scale_shape_manual(values=c(1, 2, 3, 4), labels=c("Q1", "Q2", "Q3", "Q4"))+
  stat_stars(aes(colour=Y), alpha=0.5)+
  stat_ellipse((aes(colour=Y)), level=0.95, linewidth = 1.5, type="t", alpha=0.5)
ggsave("SQePLSDA500clMen.png", width = 10, height = 10, dpi=320)

#PCA, top 500 probes closest to menarche
Y.pca <- pca(betasMentop500_Cl, ncomp=463) 

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
  labs(title= "c) PCA, top 500, closest", colour="Time to Menarche", shape="Time to Menarche")+
  xlab(label= paste("PC1 (expl. variance= ", PCA1_var, "%)",sep=""))+
  ylab(label= paste("PC2 (expl. variance= ", PCA2_var, "%)",sep=""))+
  theme(legend.position = "bottom")+
  scale_color_manual(values = c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728"), labels=c("Q1", "Q2", "Q3", "Q4"))+
  scale_shape_manual(values=c(1, 2, 3, 4), labels=c("Q1", "Q2", "Q3", "Q4"))+
  stat_stars(aes(colour=Y), alpha=0.5)+
  stat_ellipse((aes(colour=Y)), level=0.95, linewidth = 1.5, type="t", alpha=0.5)
ggsave("SQcPCA500clMen.png", width = 10, height = 10, dpi=320)

#Furthest from menarche
Top500Men_Fu <- read.csv("~/Documents/ALSPAC Methylome/Data/DMPfMentop500.csv")
Top500Men_Fu <- Top500Men_Fu[,1]

#Subset betasMen for top 500_Fu
betasMentop500_Fu <- Mennormc[,match(Top500Men_Fu,colnames(Mennormc))]
write.csv(betasMentop500_Fu,"Top500Men_Fu_norm.csv")

#PLSDA, top 500 probes furthest from menarche
plsda.dataY<- plsda(betasMentop500_Fu, Y, ncomp=463)

#Calculate variance
ev <- explained_variance(plsda.dataY$X, plsda.dataY$variates$X, ncomp=463)
Comp1v<-round(ev[1],digits = 2)*100
Comp2v<-round(ev[2],digits = 2)*100

plsdadata_plot <- data.frame(plsda.dataY$variates$X[,1:2])

ggplot(plsdadata_plot,aes(comp1,comp2,colour=Y, shape=Y))+
  geom_point(size=5, alpha=0.8)+
  coord_fixed(ratio=1)+
  theme_classic(base_size=30)+
  labs(title= "f) PLSDA, top 500, furthest", colour="Time to Menarche", shape="Time to Menarche")+
  xlab(label= paste("Component 1 (expl. variance= ", Comp1v, "%)",sep=""))+
  ylab(label= paste("Component 2 (expl. variance= ", Comp2v, "%)",sep=""))+
  theme(legend.position= "bottom")+
  scale_color_manual(values = c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728"), labels=c("Q1", "Q2", "Q3", "Q4"))+
  scale_shape_manual(values=c(1, 2, 3, 4), labels=c("Q1", "Q2", "Q3", "Q4"))+
  stat_stars(aes(colour=Y), alpha=0.5)+
  stat_ellipse((aes(colour=Y)), level=0.95, linewidth = 1.5, type="t", alpha=0.5)
ggsave("SQfPLSDA500fuMen.png", width = 10, height = 10, dpi=320)

#PCA, top 500 probes furthest from menarche
Y.pca <- pca(betasMentop500_Fu, ncomp=463) 

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
  labs(title= "d) PCA, top 500, furthest", colour="Time to Menarche", shape="Time to Menarche")+
  xlab(label= paste("PC1 (expl. variance= ", PCA1_var, "%)",sep=""))+
  ylab(label= paste("PC2 (expl. variance= ", PCA2_var, "%)",sep=""))+
  theme(legend.position = "bottom")+
  scale_color_manual(values = c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728"), labels=c("Q1", "Q2", "Q3", "Q4"))+
  scale_shape_manual(values=c(1, 2, 3, 4), labels=c("Q1", "Q2", "Q3", "Q4"))+
  stat_stars(aes(colour=Y), alpha=0.5)+
  stat_ellipse((aes(colour=Y)), level=0.95, linewidth = 1.5, type="t", alpha=0.5)
ggsave("SQdPCA500fuMen.png", width = 10, height = 10, dpi=320)

