library(mixOmics)
library(rgl)
library(tidyverse)
library(ggpubr)
library(limma)

sample_info <- read_csv("~/Documents/Code/RMCH_MethArrayBatch1/RMCHMethChAMP/Data/Meth1to46_prelim.csv")

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

#PCA all MF ##----
pcadata_all<-pca(MyNormt,ncomp=45)
Y <- as.factor(group)

#Calculate variance
eigs <- pcadata_all$sdev^2
PCA1_all_var <- eigs[1] / sum(eigs)
PCA1_all_var<-round(PCA1_all_var,digits = 2)*100

PCA2_all_var <- eigs[2] / sum(eigs)
PCA2_all_var<-round(PCA2_all_var,digits = 2)*100

pcadata_plot <- data.frame(pcadata_all$variates$X[,1:2])

ggplot(pcadata_plot,aes(PC1,PC2,colour=Y, shape=Y))+
  geom_point(size=7)+
  coord_fixed(ratio=1)+
  theme_classic(base_size=30)+
  labs(title= "a) PCA, all probes, MF", colour="Pubertal Status", shape="Pubertal Status")+
  xlab(label= paste("PC1 (expl. variance= ", PCA1_all_var, "%)",sep=""))+
  ylab(label= paste("PC2 (expl. variance= ", PCA2_all_var, "%)",sep=""))+
  theme(legend.position = "bottom")+
  scale_color_manual(values = c("#1F77B4", "#FF7F0E"))+
  stat_stars(aes(colour=Y))+
  stat_ellipse((aes(colour=Y)), level=0.95, linewidth = 1.5, type="t")
ggsave("aPCAallprobesMFpuberty_Meth.png", width = 10, height = 10, dpi=320)

Z <- as.factor(sample_info$Sex)
ggplot(pcadata_plot,aes(PC1,PC2,colour=Z, shape=Z))+
  geom_point(size=7)+
  coord_fixed(ratio=1)+
  theme_classic(base_size=30)+
  labs(title= "b) PCA, all probes, MF", colour="Sex", shape="Sex")+
  xlab(label= paste("PC1 (expl. variance= ", PCA1_all_var, "%)",sep=""))+
  ylab(label= paste("PC2 (expl. variance= ", PCA2_all_var, "%)",sep=""))+
  theme(legend.position = "bottom")+
  scale_color_manual(values = c("#2CA02C", "#D62728"))+
  stat_stars(aes(colour=Z))+
  stat_ellipse((aes(colour=Z)), level=0.95, linewidth = 1.5, type="t")
ggsave("bPCAallprobesMFsex_Meth.png", width = 10, height = 10, dpi=320)

A <- as.factor(sample_info$Ethnicity)
ggplot(pcadata_plot,aes(PC1,PC2,colour=A, shape=A))+
  geom_point(size=7)+
  coord_fixed(ratio=1)+
  theme_classic(base_size=30)+
  labs(title= "c) PCA, all probes, MF", colour="Ethnic Group", shape="Ethnic Group")+
  xlab(label= paste("PC1 (expl. variance= ", PCA1_all_var, "%)",sep=""))+
  ylab(label= paste("PC2 (expl. variance= ", PCA2_all_var, "%)",sep=""))+
  theme(legend.position = "bottom")+
  scale_color_manual(values = c("#9467BD", "#8C564B"))+
  stat_stars(aes(colour=A))+
  stat_ellipse((aes(colour=A)), level=0.95, linewidth = 1.5, type="t")
ggsave("cPCAallprobesMFethnicity_Meth.png", width = 10, height = 10, dpi=320)

B <- as.factor(sample_info$AgeGroups)
ggplot(pcadata_plot,aes(PC1,PC2,colour=B, shape=B))+
  geom_point(size=7)+
  coord_fixed(ratio=1)+
  theme_classic(base_size=30)+
  labs(title= "d) PCA, all probes, MF", colour="Age Group", shape="Age Group")+
  xlab(label= paste("PC1 (expl. variance= ", PCA1_all_var, "%)",sep=""))+
  ylab(label= paste("PC2 (expl. variance= ", PCA2_all_var, "%)",sep=""))+
  theme(legend.position = "bottom")+
  scale_color_manual(values = c("#E377C2", "#7F7F7F"))+
  stat_stars(aes(colour=B))+
  stat_ellipse((aes(colour=B)), level=0.95, linewidth = 1.5, type="t")
ggsave("dPCAallprobesMFagegp_Meth.png", width = 10, height = 10, dpi=320)

#Explained variance
plot(pcadata_all)
plotLoadings(pcadata_all)


#PLSDA all MF----
plsda_all<- plsda(MyNormt, Y, ncomp = 45)

#Calculate variance
ev <- explained_variance(plsda_all$X, plsda_all$variates$X, ncomp =45)
Comp1v<-round(ev[1],digits = 2)*100
Comp2v<-round(ev[2],digits = 2)*100

plsdadata_plot <- data.frame(plsda_all$variates$X[,1:2])

ggplot(plsdadata_plot,aes(comp1,comp2,colour=Y, shape=Y))+
  geom_point(size=7)+
  coord_fixed(ratio = 1)+
  theme_classic(base_size=30)+
  labs(title= "e) PLSDA, all probes, MF", colour="Pubertal Status", shape="Pubertal Status")+
  xlab(label= paste("Component 1 (expl. variance= ", Comp1v, "%)",sep=""))+
  ylab(label= paste("Component 2 (expl. variance= ", Comp2v, "%)",sep=""))+
  theme(legend.position= "bottom")+
  scale_color_manual(values = c("#1F77B4", "#FF7F0E"))+
  stat_stars(aes(colour=Y))+
  stat_ellipse((aes(colour=Y)), level=0.95, linewidth = 1.5, type="t")
ggsave("ePLSDAallgenesMFpuberty_Meth.png", width = 10, height = 10, dpi=320)

plsda_all<- plsda(MyNormt, Z, ncomp = 45)
ev <- explained_variance(plsda_all$X, plsda_all$variates$X, ncomp =45)
Comp1v<-round(ev[1],digits = 2)*100
Comp2v<-round(ev[2],digits = 2)*100
plsdadata_plot <- data.frame(plsda_all$variates$X[,1:2])

ggplot(plsdadata_plot,aes(comp1,comp2,colour=Z, shape=Z))+
  geom_point(size=7)+
  coord_fixed(ratio = 1)+
  theme_classic(base_size=30)+
  labs(title= "f) PLSDA, all probes, MF", colour="Sex", shape="Sex")+
  xlab(label= paste("Component 1 (expl. variance= ", Comp1v, "%)",sep=""))+
  ylab(label= paste("Component 2 (expl. variance= ", Comp2v, "%)",sep=""))+
  theme(legend.position= "bottom")+
  scale_color_manual(values = c("#2CA02C", "#D62728"))+
  stat_stars(aes(colour=Z))+
  stat_ellipse((aes(colour=Z)), level=0.95, linewidth = 1.5, type="t")
ggsave("fPLSDAallprobesMFsex_Meth.png", width = 10, height = 10, dpi=320)

plsda_all<- plsda(MyNormt, A, ncomp = 45)
ev <- explained_variance(plsda_all$X, plsda_all$variates$X, ncomp =45)
Comp1v<-round(ev[1],digits = 2)*100
Comp2v<-round(ev[2],digits = 2)*100
plsdadata_plot <- data.frame(plsda_all$variates$X[,1:2])

ggplot(plsdadata_plot,aes(comp1,comp2,colour=A, shape=A))+
  geom_point(size=7)+
  coord_fixed(ratio = 1)+
  theme_classic(base_size=30)+
  labs(title= "g) PLSDA, all probes, MF", colour="Ethnic Group", shape="Ethnic Group")+
  xlab(label= paste("Component 1 (expl. variance= ", Comp1v, "%)",sep=""))+
  ylab(label= paste("Component 2 (expl. variance= ", Comp2v, "%)",sep=""))+
  theme(legend.position= "bottom")+
  scale_color_manual(values = c("#9467BD", "#8C564B"))+
  stat_stars(aes(colour=A))+
  stat_ellipse((aes(colour=A)), level=0.95, linewidth = 1.5, type="t")
ggsave("gPLSDAallprobesMFethnicity_Meth.png", width = 10, height = 10, dpi=320)


plsda_all<- plsda(MyNormt, B, ncomp = 45)
ev <- explained_variance(plsda_all$X, plsda_all$variates$X, ncomp =45)
Comp1v<-round(ev[1],digits = 2)*100
Comp2v<-round(ev[2],digits = 2)*100
plsdadata_plot <- data.frame(plsda_all$variates$X[,1:2])

ggplot(plsdadata_plot,aes(comp1,comp2,colour=B, shape=B))+
  geom_point(size=7)+
  coord_fixed(ratio = 1)+
  theme_classic(base_size=30)+
  labs(title= "h) PLSDA, all probes, MF", colour="Age Group", shape="Age Group")+
  xlab(label= paste("Component 1 (expl. variance= ", Comp1v, "%)",sep=""))+
  ylab(label= paste("Component 2 (expl. variance= ", Comp2v, "%)",sep=""))+
  theme(legend.position= "bottom")+
  scale_color_manual(values = c("#E377C2", "#7F7F7F"))+
  stat_stars(aes(colour=B))+
  stat_ellipse((aes(colour=B)), level=0.95, linewidth = 1.5, type="t")
ggsave("hPLSDAallprobesMFagegp_Meth.png", width = 10, height = 10, dpi=320)

#Create subset of normalised betas for top 500 differential probes
Top500 <- rownames(DMPtop500)
Top500normbeta<- MyNormt[,match(Top500,colnames(MyNormt))]
  
#PCA top 500 - contains normalised counts for all the data MF----
pcadata_500<- pca(Top500normbeta, ncomp=45)
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
  labs(title= "a) PCA, top 500, MF", colour="Pubertal Status", shape="Pubertal Status")+
  xlab(label= paste("PC1 (expl. variance= ", PCA1_500_var, "%)",sep=""))+
  ylab(label= paste("PC2 (expl. variance= ", PCA2_500_var, "%)",sep=""))+
  theme(legend.position = "bottom")+
  scale_color_manual(values = c("#1F77B4", "#FF7F0E"))+
  stat_stars(aes(colour=Y))+
  stat_ellipse((aes(colour=Y)), level=0.95, linewidth = 1.5, type="t")
ggsave("aPCAtop500MFpuberty_Meth.png", width = 10, height = 10, dpi=320)

#Explained variance
plot(pcadata_500)
plotLoadings(pcadata_500)

#PLSDA top 500 MF----
plsda_500<- plsda(Top500normbeta, Y, ncomp = 45)

#Calculate variance
ev <- explained_variance(plsda_500$X, plsda_500$variates$X, ncomp =45)
Comp1v<-round(ev[1],digits = 2)*100
Comp2v<-round(ev[2],digits = 2)*100

plsdadata_plot <- data.frame(plsda_500$variates$X[,1:2])

ggplot(plsdadata_plot,aes(comp1,comp2,colour=Y, shape=Y))+
  geom_point(size=7)+
  coord_fixed(ratio = 1)+
  theme_classic(base_size=30)+
  labs(title= "b) PLSDA, top 500, MF", colour="Pubertal Status", shape="Pubertal Status")+
  xlab(label= paste("Component 1 (expl. variance= ", Comp1v, "%)",sep=""))+
  ylab(label= paste("Component 2 (expl. variance= ", Comp2v, "%)",sep=""))+
  theme(legend.position= "bottom")+
  scale_color_manual(values = c("#1F77B4", "#FF7F0E"))+
  stat_stars(aes(colour=Y))+
  stat_ellipse((aes(colour=Y)), level=0.95, linewidth = 1.5, type="t")
ggsave("bPLSDAtop500MFpuberty_Meth.png", width = 10, height = 10, dpi=320)

#RF MF
library(randomForest)
library(Boruta)
library(pROC)
library(ROCR)
library(smotefamily)
library(caret)
library(mltools)
library(ggrepel)

#Divide into pubertal and pre-pubertal
clindata_prepub <- sample_info[sample_info$Group== "Prepubertal",]
clindata_pub <- sample_info[sample_info$Group== "Pubertal",]

#Randomly Split Dataset
set.seed(123) 
n <- nrow(clindata_prepub)  # Number of rows in the dataframe
split_ids <- sample(rep(1:4, length.out = n))  # Assign random numbers 1 to 4
clindata_prepub_parts <- split(clindata_prepub, split_ids)  # Split dataframe into 4 parts

n <- nrow(clindata_pub)  # Number of rows in the dataframe
split_ids <- sample(rep(1:4, length.out = n))  # Assign random numbers 1 to 4
clindata_pub_parts <- split(clindata_pub, split_ids)  # Split dataframe into 4 parts

Part1 <- rbind(clindata_prepub_parts[["1"]],clindata_pub_parts[["4"]])
Part2 <- rbind(clindata_prepub_parts[["2"]],clindata_pub_parts[["2"]])
Part3 <- rbind(clindata_prepub_parts[["3"]],clindata_pub_parts[["3"]])
Part4 <- rbind(clindata_prepub_parts[["4"]],clindata_pub_parts[["1"]])

Comb1 <- rbind(Part1,Part2,Part3)
Comb2 <- rbind(Part1,Part2,Part4)
Comb3 <- rbind(Part1,Part3,Part4)
Comb4 <- rbind(Part2,Part3,Part4)

Top500_Part1 <- as.data.frame(Top500normbeta[match(Part1$ParticipantID,rownames(Top500normbeta)),])
Top500_Part2 <- as.data.frame(Top500normbeta[match(Part2$ParticipantID,rownames(Top500normbeta)),])
Top500_Part3 <- as.data.frame(Top500normbeta[match(Part3$ParticipantID,rownames(Top500normbeta)),])
Top500_Part4 <- as.data.frame(Top500normbeta[match(Part4$ParticipantID,rownames(Top500normbeta)),])
Top500_Comb1 <- as.data.frame(Top500normbeta[match(Comb1$ParticipantID,rownames(Top500normbeta)),])
Top500_Comb2 <- as.data.frame(Top500normbeta[match(Comb2$ParticipantID,rownames(Top500normbeta)),])
Top500_Comb3 <- as.data.frame(Top500normbeta[match(Comb3$ParticipantID,rownames(Top500normbeta)),])
Top500_Comb4 <- as.data.frame(Top500normbeta[match(Comb4$ParticipantID,rownames(Top500normbeta)),])

#Boruta & RF - Comb1 (Parts 1,2,3): train, Part 4: test----
set.seed (123)
pred.var <- Comb1$Group
data.train<-cbind(Top500_Comb1,pred.var)

#Boruta
bor<-Boruta(pred.var~.,data=data.train,maxRuns=1000)
bor
plot(bor)
Boruta_Comb1_imphist <- bor[["ImpHistory"]]
write.csv(bor[["ImpHistory"]],file="Boruta_Comb1_imphistoryfull.csv")
write.csv(bor[["finalDecision"]],file="Boruta_Comb1.csv")
imp.probes_Comb1<-names(bor$finalDecision[which(bor$finalDecision=="Confirmed")])

data.train<-data.frame(cbind(Top500_Comb1[match(imp.probes_Comb1,colnames(Top500_Comb1))],pred.var))

#Running RF
rf<-randomForest(pred.var~.,data=data.train,ntree=1000)
rf

#OOB AUC
OOBvotes_Comb1=predict(rf,type = "prob")
write.csv((data.frame(OOBvotes_Comb1)),file="OOBvotes_Comb1_full.csv")

OOBpred=prediction(OOBvotes_Comb1[,2],data.train$pred.var)
# Area under curve
OOBauc=performance(OOBpred, "auc")
OOBauc <- OOBauc@y.values[[1]]
OOBauc
# True Positive and Negative Rate
OOBtpfp=performance(OOBpred, "tpr","fpr")

#Preparation of test data (unsupervised)
pred.var<-as.factor(Part4$Group) #predictive variable
data.test<- Top500_Part4[,match(imp.probes_Comb1,colnames(Top500_Part4))]

#Validation cohort 
Valpredgps<-predict(rf,data.test)
ConMat <- confusionMatrix(Valpredgps,pred.var)
ConMat
F1 <- ConMat[["byClass"]][["F1"]]
F1
Table <- as.vector(ConMat[["table"]])
Tablem <-matrix(c(Table),nrow=2)
MatCC <- mcc(confusionM=Tablem)
MatCC

#Validation AUC
Valvotes_Comb1=predict(rf,newdata = data.test,type = "prob")
write.csv((data.frame(Valvotes_Comb1)),file="Valvotes_Comb1.csv")
Valpred=prediction(Valvotes_Comb1[,2],pred.var)
#Area under curve
Valauc=performance(Valpred, "auc")
Valauc <- Valauc@y.values[[1]]
Valauc 
#True Positive and Negative Rate
Valtpfp=performance(Valpred, "tpr","fpr")

#Results
write.csv((data.frame(OOBauc,Valauc,F1, MatCC)),file="Perf_Comb1.csv")

OFPR <- unlist(OOBtpfp@x.values)
OTPR <- unlist(OOBtpfp@y.values)
OROC <- data.frame(OFPR,OTPR)
write.csv(OROC,file="OROC_Comb1.csv")

VFPR <- unlist(Valtpfp@x.values)
VTPR <- unlist(Valtpfp@y.values)
VROC <- data.frame(VFPR,VTPR)
write.csv(VROC,file="VROC_Comb1.csv")

#Boruta & RF - Comb2 (Parts 1,2,4): train, Part 3: test----
set.seed (123)
pred.var <- Comb2$Group
data.train<-cbind(Top500_Comb2,pred.var)

#Boruta
bor<-Boruta(pred.var~.,data=data.train,maxRuns=1000)
bor
plot(bor)
Boruta_Comb2_imphist <- bor[["ImpHistory"]]
write.csv(bor[["ImpHistory"]],file="Boruta_Comb2_imphistoryfull.csv")
write.csv(bor[["finalDecision"]],file="Boruta_Comb2.csv")
imp.probes_Comb2<-names(bor$finalDecision[which(bor$finalDecision=="Confirmed")])

data.train<-data.frame(cbind(Top500_Comb2[match(imp.probes_Comb2,colnames(Top500_Comb2))],pred.var))

#Running RF
rf<-randomForest(pred.var~.,data=data.train,ntree=1000)
rf

#OOB AUC
OOBvotes_Comb2=predict(rf,type = "prob")
write.csv((data.frame(OOBvotes_Comb2)),file="OOBvotes_Comb2.csv")
OOBpred=prediction(OOBvotes_Comb2[,2],data.train$pred.var)
# Area under curve
OOBauc=performance(OOBpred, "auc")
OOBauc <- OOBauc@y.values[[1]]
OOBauc
# True Positive and Negative Rate
OOBtpfp=performance(OOBpred, "tpr","fpr")

#Preparation of test data (unsupervised)
pred.var<-as.factor(Part3$Group) #predictive variable
data.test<- Top500_Part3[,match(imp.probes_Comb2,colnames(Top500_Part3))]

#Validation cohort 
Valpredgps<-predict(rf,data.test)
ConMat <- confusionMatrix(Valpredgps,pred.var)
ConMat
F1 <- ConMat[["byClass"]][["F1"]]
F1
Table <- as.vector(ConMat[["table"]])
Tablem <-matrix(c(Table),nrow=2)
MatCC <- mcc(confusionM=Tablem)
MatCC

#Validation AUC
Valvotes_Comb2=predict(rf,newdata = data.test,type = "prob")
write.csv((data.frame(Valvotes_Comb2)),file="Valvotes_Comb2.csv")
Valpred=prediction(Valvotes_Comb2[,2],pred.var)
#Area under curve
Valauc=performance(Valpred, "auc")
Valauc <- Valauc@y.values[[1]]
Valauc 
#True Positive and Negative Rate
Valtpfp=performance(Valpred, "tpr","fpr")

#Results
write.csv((data.frame(OOBauc,Valauc,F1, MatCC)),file="Perf_Comb2.csv")

OFPR <- unlist(OOBtpfp@x.values)
OTPR <- unlist(OOBtpfp@y.values)
OROC <- data.frame(OFPR,OTPR)
write.csv(OROC,file="OROC_Comb2.csv")

VFPR <- unlist(Valtpfp@x.values)
VTPR <- unlist(Valtpfp@y.values)
VROC <- data.frame(VFPR,VTPR)
write.csv(VROC,file="VROC_Comb2.csv")

summary(data.train$pred.var)
summary(pred.var)

#Boruta & RF - Comb3 (Parts 1,3,4): train, Part 2: test----
set.seed (123)
pred.var <- Comb3$Group
data.train<-cbind(Top500_Comb3,pred.var)

#Boruta
bor<-Boruta(pred.var~.,data=data.train,maxRuns=1000)
bor
plot(bor)
Boruta_Comb3_imphist <- bor[["ImpHistory"]]
write.csv(bor[["ImpHistory"]],file="Boruta_Comb3_imphistoryfull.csv")
write.csv(bor[["finalDecision"]],file="Boruta_Comb3.csv")
imp.probes_Comb3<-names(bor$finalDecision[which(bor$finalDecision=="Confirmed")])

data.train<-data.frame(cbind(Top500_Comb3[match(imp.probes_Comb3,colnames(Top500_Comb3))],pred.var))

#Running RF
rf<-randomForest(pred.var~.,data=data.train,ntree=1000)
rf

#OOB AUC
OOBvotes_Comb3=predict(rf,type = "prob")
write.csv((data.frame(OOBvotes_Comb3)),file="OOBvotes_Comb3.csv")
OOBpred=prediction(OOBvotes_Comb3[,2],data.train$pred.var)
# Area under curve
OOBauc=performance(OOBpred, "auc")
OOBauc <- OOBauc@y.values[[1]]
OOBauc
# True Positive and Negative Rate
OOBtpfp=performance(OOBpred, "tpr","fpr")

#Preparation of test data (unsupervised)
pred.var<-as.factor(Part2$Group) #predictive variable
data.test<- Top500_Part2[,match(imp.probes_Comb3,colnames(Top500_Part2))]

#Validation cohort 
Valpredgps<-predict(rf,data.test)
ConMat <- confusionMatrix(Valpredgps,pred.var)
ConMat
F1 <- ConMat[["byClass"]][["F1"]]
F1
Table <- as.vector(ConMat[["table"]])
Tablem <-matrix(c(Table),nrow=2)
MatCC <- mcc(confusionM=Tablem)
MatCC

#Validation AUC
Valvotes_Comb3=predict(rf,newdata = data.test,type = "prob")
write.csv((data.frame(Valvotes_Comb3)),file="Valvotes_Comb3.csv")
Valpred=prediction(Valvotes_Comb3[,2],pred.var)
#Area under curve
Valauc=performance(Valpred, "auc")
Valauc <- Valauc@y.values[[1]]
Valauc 
#True Positive and Negative Rate
Valtpfp=performance(Valpred, "tpr","fpr")

#Results
write.csv((data.frame(OOBauc,Valauc,F1, MatCC)),file="Perf_Comb3.csv")

OFPR <- unlist(OOBtpfp@x.values)
OTPR <- unlist(OOBtpfp@y.values)
OROC <- data.frame(OFPR,OTPR)
write.csv(OROC,file="OROC_Comb3.csv")

VFPR <- unlist(Valtpfp@x.values)
VTPR <- unlist(Valtpfp@y.values)
VROC <- data.frame(VFPR,VTPR)
write.csv(VROC,file="VROC_Comb3.csv")

#Boruta & RF - Comb4 (Parts 2,3,4): train, Part 1: test----
set.seed (123)
pred.var <- Comb4$Group
data.train<-cbind(Top500_Comb4,pred.var)

#Boruta
bor<-Boruta(pred.var~.,data=data.train,maxRuns=1000)
bor
plot(bor)
Boruta_Comb4_imphist <- bor[["ImpHistory"]]
write.csv(bor[["ImpHistory"]],file="Boruta_Comb4_imphistoryfull.csv")
write.csv(bor[["finalDecision"]],file="Boruta_Comb4.csv")
imp.probes_Comb4<-names(bor$finalDecision[which(bor$finalDecision=="Confirmed")])

data.train<-data.frame(cbind(Top500_Comb4[match(imp.probes_Comb4,colnames(Top500_Comb4))],pred.var))

#Running RF
rf<-randomForest(pred.var~.,data=data.train,ntree=1000)
rf

#OOB AUC
OOBvotes_Comb4=predict(rf,type = "prob")
write.csv((data.frame(OOBvotes_Comb4)),file="OOBvotes_Comb4.csv")
OOBpred=prediction(OOBvotes_Comb4[,2],data.train$pred.var)
# Area under curve
OOBauc=performance(OOBpred, "auc")
OOBauc <- OOBauc@y.values[[1]]
OOBauc
# True Positive and Negative Rate
OOBtpfp=performance(OOBpred, "tpr","fpr")

#Preparation of test data (unsupervised)
pred.var<-as.factor(Part1$Group) #predictive variable
data.test<- Top500_Part1[,match(imp.probes_Comb4,colnames(Top500_Part1))]

#Validation cohort 
Valpredgps<-predict(rf,data.test)
ConMat <- confusionMatrix(Valpredgps,pred.var)
ConMat
F1 <- ConMat[["byClass"]][["F1"]]
F1
Table <- as.vector(ConMat[["table"]])
Tablem <-matrix(c(Table),nrow=2)
MatCC <- mcc(confusionM=Tablem)
MatCC

#Validation AUC
Valvotes_Comb4=predict(rf,newdata = data.test,type = "prob")
write.csv((data.frame(Valvotes_Comb4)),file="Valvotes_Comb4.csv")
Valpred=prediction(Valvotes_Comb4[,2],pred.var)
#Area under curve
Valauc=performance(Valpred, "auc")
Valauc <- Valauc@y.values[[1]]
Valauc 
#True Positive and Negative Rate
Valtpfp=performance(Valpred, "tpr","fpr")

#Results
write.csv((data.frame(OOBauc,Valauc,F1, MatCC)),file="Perf_Comb4.csv")

OFPR <- unlist(OOBtpfp@x.values)
OTPR <- unlist(OOBtpfp@y.values)
OROC <- data.frame(OFPR,OTPR)
write.csv(OROC,file="OROC_Comb4.csv")

VFPR <- unlist(Valtpfp@x.values)
VTPR <- unlist(Valtpfp@y.values)
VROC <- data.frame(VFPR,VTPR)
write.csv(VROC,file="VROC_Comb4.csv")

#Individual Box plots (MF)----
topprobes_betas <- Top500normbeta[,match(topprobes,colnames(Top500normbeta))]
topprobes_betas <- as.data.frame(topprobes_betas) 
topprobes_betas$SampleID <- rownames(topprobes_betas)
topprobes_betas <- merge(topprobes_betas, sample_info, by = "SampleID", all = TRUE)

ggplot(topprobes_betas, aes(x =Group, y =`cg13373080_BC21`, fill=Group)) +
  geom_boxplot()+
  theme_classic(base_size=14)+
  theme(legend.position ="none")+
  labs(x = NULL, y = "cg13373080_BC21 Methylation (normalised beta values)")
ggsave("cg13373080_BC21.png", height=5, width=2.5, dpi=320)

#etc.

#DMP info for predictive probes----
DMP_pred <- DMPtop500[match(imp.probes_all,rownames(DMPtop500)),]
write.csv(DMP_pred, file="DMP_Pred.csv")

imp.probes_incdups <- c(imp.probes_Comb1,imp.probes_Comb2, imp.probes_Comb3, imp.probes_Comb4)
DMP_pred_incdups <- DMPtop500[match(imp.probes_incdups,rownames(DMPtop500)),]
write.csv(DMP_pred_incdups, file="DMP_Pred_incdups.csv")

#Females only----
summary(sample_info$Sex)
sample_infoF <- sample_info[sample_info$Sex=="Female",]

#Create comparison groups----
group <- sample_infoF$Group
ethnicity <- sample_infoF$Ethnicity
agegp <- sample_infoF$AgeGroups
design <- model.matrix(~0+group+ethnicity+agegp)
my_contrast <- makeContrasts(groupPubertal - groupPrepubertal,levels= design)

#Load files
load("~/Documents/Code/RMCH_MethArrayBatch1/RMCHMethChAMP/MyNormFt.Robj")

#PCA all F ##----
pcadata_all<-pca(MyNormFt,ncomp=31)
Y <- as.factor(group)

#Calculate variance
eigs <- pcadata_all$sdev^2
PCA1_all_var <- eigs[1] / sum(eigs)
PCA1_all_var<-round(PCA1_all_var,digits = 2)*100

PCA2_all_var <- eigs[2] / sum(eigs)
PCA2_all_var<-round(PCA2_all_var,digits = 2)*100

pcadata_plot <- data.frame(pcadata_all$variates$X[,1:2])

ggplot(pcadata_plot,aes(PC1,PC2,colour=Y, shape=Y))+
  geom_point(size=7)+
  coord_fixed(ratio=1)+
  theme_classic(base_size=30)+
  labs(title= "a) PCA, all probes, F", colour="Pubertal Status", shape="Pubertal Status")+
  xlab(label= paste("PC1 (expl. variance= ", PCA1_all_var, "%)",sep=""))+
  ylab(label= paste("PC2 (expl. variance= ", PCA2_all_var, "%)",sep=""))+
  theme(legend.position = "bottom")+
  scale_color_manual(values = c("#1F77B4", "#FF7F0E"))+
  stat_stars(aes(colour=Y))+
  stat_ellipse((aes(colour=Y)), level=0.95, linewidth = 1.5, type="t")
ggsave("aPCAallprobesFpuberty_Meth.png", width = 10, height = 10, dpi=320)

A <- as.factor(sample_infoF$Ethnicity)
ggplot(pcadata_plot,aes(PC1,PC2,colour=A, shape=A))+
  geom_point(size=7)+
  coord_fixed(ratio=1)+
  theme_classic(base_size=30)+
  labs(title= "b) PCA, all probes, F", colour="Ethnic Group", shape="Ethnic Group")+
  xlab(label= paste("PC1 (expl. variance= ", PCA1_all_var, "%)",sep=""))+
  ylab(label= paste("PC2 (expl. variance= ", PCA2_all_var, "%)",sep=""))+
  theme(legend.position = "bottom")+
  scale_color_manual(values = c("#9467BD", "#8C564B"))+
  stat_stars(aes(colour=A))+
  stat_ellipse((aes(colour=A)), level=0.95, linewidth = 1.5, type="t")
ggsave("bPCAallprobesFethnicity_Meth.png", width = 10, height = 10, dpi=320)

#Explained variance
plot(pcadata_all)
plotLoadings(pcadata_all)

#PLSDA all F----
plsda_all<- plsda(MyNormFt, Y, ncomp = 31)

#Calculate variance
ev <- explained_variance(plsda_all$X, plsda_all$variates$X, ncomp =31)
Comp1v<-round(ev[1],digits = 2)*100
Comp2v<-round(ev[2],digits = 2)*100

plsdadata_plot <- data.frame(plsda_all$variates$X[,1:2])

ggplot(plsdadata_plot,aes(comp1,comp2,colour=Y, shape=Y))+
  geom_point(size=7)+
  coord_fixed(ratio = 1)+
  theme_classic(base_size=30)+
  labs(title= "c) PLSDA, all probes, F", colour="Pubertal Status", shape="Pubertal Status")+
  xlab(label= paste("Component 1 (expl. variance= ", Comp1v, "%)",sep=""))+
  ylab(label= paste("Component 2 (expl. variance= ", Comp2v, "%)",sep=""))+
  theme(legend.position= "bottom")+
  scale_color_manual(values = c("#1F77B4", "#FF7F0E"))+
  stat_stars(aes(colour=Y))+
  stat_ellipse((aes(colour=Y)), level=0.95, linewidth = 1.5, type="t")
ggsave("cPLSDAallgenesFpuberty_Meth.png", width = 10, height = 10, dpi=320)

plsda_all<- plsda(MyNormFt, A, ncomp = 31)
ev <- explained_variance(plsda_all$X, plsda_all$variates$X, ncomp =31)
Comp1v<-round(ev[1],digits = 2)*100
Comp2v<-round(ev[2],digits = 2)*100
plsdadata_plot <- data.frame(plsda_all$variates$X[,1:2])

ggplot(plsdadata_plot,aes(comp1,comp2,colour=A, shape=A))+
  geom_point(size=7)+
  coord_fixed(ratio = 1)+
  theme_classic(base_size=30)+
  labs(title= "d) PLSDA, all probes, F", colour="Ethnic Group", shape="Ethnic Group")+
  xlab(label= paste("Component 1 (expl. variance= ", Comp1v, "%)",sep=""))+
  ylab(label= paste("Component 2 (expl. variance= ", Comp2v, "%)",sep=""))+
  theme(legend.position= "bottom")+
  scale_color_manual(values = c("#9467BD", "#8C564B"))+
  stat_stars(aes(colour=A))+
  stat_ellipse((aes(colour=A)), level=0.95, linewidth = 1.5, type="t")
ggsave("dPLSDAallprobesFethnicity_Meth.png", width = 10, height = 10, dpi=320)

#Create subset of normalised betas for top 500 differential probes
DMPtop500F <- read.csv("~/Documents/Code/RMCH_MethArrayBatch1/RMCHMethChAMP/DMPtop500F.csv")
Top500 <- DMPtop500F$X
Top500normbeta<- MyNormFt[,match(Top500,colnames(MyNormFt))]

#PCA top 500 - contains normalised counts for all the data F----
pcadata_500<- pca(Top500normbeta, ncomp=31)
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
  labs(title= "c) PCA, top 500, F", colour="Pubertal Status", shape="Pubertal Status")+
  xlab(label= paste("PC1 (expl. variance= ", PCA1_500_var, "%)",sep=""))+
  ylab(label= paste("PC2 (expl. variance= ", PCA2_500_var, "%)",sep=""))+
  theme(legend.position = "bottom")+
  scale_color_manual(values = c("#1F77B4", "#FF7F0E"))+
  stat_stars(aes(colour=Y))+
  stat_ellipse((aes(colour=Y)), level=0.95, linewidth = 1.5, type="t")
ggsave("cPCAtop500Fpuberty_Meth.png", width = 10, height = 10, dpi=320)

#Explained variance
plot(pcadata_500)
plotLoadings(pcadata_500)

#PLSDA top 500 MF----
plsda_500<- plsda(Top500normbeta, Y, ncomp = 31)

#Calculate variance
ev <- explained_variance(plsda_500$X, plsda_500$variates$X, ncomp =31)
Comp1v<-round(ev[1],digits = 2)*100
Comp2v<-round(ev[2],digits = 2)*100

plsdadata_plot <- data.frame(plsda_500$variates$X[,1:2])

ggplot(plsdadata_plot,aes(comp1,comp2,colour=Y, shape=Y))+
  geom_point(size=7)+
  coord_fixed(ratio = 1)+
  theme_classic(base_size=30)+
  labs(title= "d) PLSDA, top 500, F", colour="Pubertal Status", shape="Pubertal Status")+
  xlab(label= paste("Component 1 (expl. variance= ", Comp1v, "%)",sep=""))+
  ylab(label= paste("Component 2 (expl. variance= ", Comp2v, "%)",sep=""))+
  theme(legend.position= "bottom")+
  scale_color_manual(values = c("#1F77B4", "#FF7F0E"))+
  stat_stars(aes(colour=Y))+
  stat_ellipse((aes(colour=Y)), level=0.95, linewidth = 1.5, type="t")
ggsave("dPLSDAtop500Fpuberty_Meth.png", width = 10, height = 10, dpi=320)

#RF MF
library(randomForest)
library(Boruta)
library(pROC)
library(ROCR)
library(smotefamily)
library(caret)
library(mltools)
library(ggrepel)

#Divide into pubertal and pre-pubertal
clindata_prepub <- sample_infoF[sample_infoF$Group== "Prepubertal",]
clindata_pub <- sample_infoF[sample_infoF$Group== "Pubertal",]

#Randomly Split Dataset
set.seed(123) 
n <- nrow(clindata_prepub)  # Number of rows in the dataframe
split_ids <- sample(rep(1:4, length.out = n))  # Assign random numbers 1 to 4
clindata_prepub_parts <- split(clindata_prepub, split_ids)  # Split dataframe into 4 parts

n <- nrow(clindata_pub)  # Number of rows in the dataframe
split_ids <- sample(rep(1:4, length.out = n))  # Assign random numbers 1 to 4
clindata_pub_parts <- split(clindata_pub, split_ids)  # Split dataframe into 4 parts

Part1 <- rbind(clindata_prepub_parts[["1"]],clindata_pub_parts[["4"]])
Part2 <- rbind(clindata_prepub_parts[["2"]],clindata_pub_parts[["2"]])
Part3 <- rbind(clindata_prepub_parts[["3"]],clindata_pub_parts[["3"]])
Part4 <- rbind(clindata_prepub_parts[["4"]],clindata_pub_parts[["1"]])

Comb1 <- rbind(Part1,Part2,Part3)
Comb2 <- rbind(Part1,Part2,Part4)
Comb3 <- rbind(Part1,Part3,Part4)
Comb4 <- rbind(Part2,Part3,Part4)

Top500_Part1 <- as.data.frame(Top500normbeta[match(Part1$ParticipantID,rownames(Top500normbeta)),])
Top500_Part2 <- as.data.frame(Top500normbeta[match(Part2$ParticipantID,rownames(Top500normbeta)),])
Top500_Part3 <- as.data.frame(Top500normbeta[match(Part3$ParticipantID,rownames(Top500normbeta)),])
Top500_Part4 <- as.data.frame(Top500normbeta[match(Part4$ParticipantID,rownames(Top500normbeta)),])
Top500_Comb1 <- as.data.frame(Top500normbeta[match(Comb1$ParticipantID,rownames(Top500normbeta)),])
Top500_Comb2 <- as.data.frame(Top500normbeta[match(Comb2$ParticipantID,rownames(Top500normbeta)),])
Top500_Comb3 <- as.data.frame(Top500normbeta[match(Comb3$ParticipantID,rownames(Top500normbeta)),])
Top500_Comb4 <- as.data.frame(Top500normbeta[match(Comb4$ParticipantID,rownames(Top500normbeta)),])

#Boruta & RF - Comb1 (Parts 1,2,3): train, Part 4: test----
set.seed (123)
pred.var <- Comb1$Group
data.train<-cbind(Top500_Comb1,pred.var)

#Boruta
bor<-Boruta(pred.var~.,data=data.train,maxRuns=1000)
bor
plot(bor)
Boruta_Comb1_imphist <- bor[["ImpHistory"]]
write.csv(bor[["ImpHistory"]],file="Boruta_Comb1_imphistoryfullF.csv")
write.csv(bor[["finalDecision"]],file="Boruta_Comb1F.csv")
imp.probes_Comb1<-names(bor$finalDecision[which(bor$finalDecision=="Confirmed")])

data.train<-data.frame(cbind(Top500_Comb1[match(imp.probes_Comb1,colnames(Top500_Comb1))],pred.var))

#Running RF
rf<-randomForest(pred.var~.,data=data.train,ntree=1000)
rf

#OOB AUC
OOBvotes_Comb1=predict(rf,type = "prob")
write.csv((data.frame(OOBvotes_Comb1)),file="OOBvotes_Comb1_fullF.csv")

OOBpred=prediction(OOBvotes_Comb1[,2],data.train$pred.var)
# Area under curve
OOBauc=performance(OOBpred, "auc")
OOBauc <- OOBauc@y.values[[1]]
OOBauc
# True Positive and Negative Rate
OOBtpfp=performance(OOBpred, "tpr","fpr")

#Preparation of test data (unsupervised)
pred.var<-as.factor(Part4$Group) #predictive variable
data.test<- Top500_Part4[,match(imp.probes_Comb1,colnames(Top500_Part4))]

#Validation cohort 
Valpredgps<-predict(rf,data.test)
ConMat <- confusionMatrix(Valpredgps,pred.var)
ConMat
F1 <- ConMat[["byClass"]][["F1"]]
F1
Table <- as.vector(ConMat[["table"]])
Tablem <-matrix(c(Table),nrow=2)
MatCC <- mcc(confusionM=Tablem)
MatCC

#Validation AUC
Valvotes_Comb1=predict(rf,newdata = data.test,type = "prob")
write.csv((data.frame(Valvotes_Comb1)),file="Valvotes_Comb1F.csv")
Valpred=prediction(Valvotes_Comb1[,2],pred.var)
#Area under curve
Valauc=performance(Valpred, "auc")
Valauc <- Valauc@y.values[[1]]
Valauc 
#True Positive and Negative Rate
Valtpfp=performance(Valpred, "tpr","fpr")

#Results
write.csv((data.frame(OOBauc,Valauc,F1, MatCC)),file="Perf_Comb1F.csv")

OFPR <- unlist(OOBtpfp@x.values)
OTPR <- unlist(OOBtpfp@y.values)
OROC <- data.frame(OFPR,OTPR)
write.csv(OROC,file="OROC_Comb1F.csv")

VFPR <- unlist(Valtpfp@x.values)
VTPR <- unlist(Valtpfp@y.values)
VROC <- data.frame(VFPR,VTPR)
write.csv(VROC,file="VROC_Comb1F.csv")

#Boruta & RF - Comb2 (Parts 1,2,4): train, Part 3: test----
set.seed (123)
pred.var <- Comb2$Group
data.train<-cbind(Top500_Comb2,pred.var)

#Boruta
bor<-Boruta(pred.var~.,data=data.train,maxRuns=1000)
bor
plot(bor)
Boruta_Comb2_imphist <- bor[["ImpHistory"]]
write.csv(bor[["ImpHistory"]],file="Boruta_Comb2_imphistoryfullF.csv")
write.csv(bor[["finalDecision"]],file="Boruta_Comb2F.csv")
imp.probes_Comb2<-names(bor$finalDecision[which(bor$finalDecision=="Confirmed")])

data.train<-data.frame(cbind(Top500_Comb2[match(imp.probes_Comb2,colnames(Top500_Comb2))],pred.var))

#Running RF
rf<-randomForest(pred.var~.,data=data.train,ntree=1000)
rf

#OOB AUC
OOBvotes_Comb2=predict(rf,type = "prob")
write.csv((data.frame(OOBvotes_Comb2)),file="OOBvotes_Comb2F.csv")
OOBpred=prediction(OOBvotes_Comb2[,2],data.train$pred.var)
# Area under curve
OOBauc=performance(OOBpred, "auc")
OOBauc <- OOBauc@y.values[[1]]
OOBauc
# True Positive and Negative Rate
OOBtpfp=performance(OOBpred, "tpr","fpr")

#Preparation of test data (unsupervised)
pred.var<-as.factor(Part3$Group) #predictive variable
data.test<- Top500_Part3[,match(imp.probes_Comb2,colnames(Top500_Part3))]

#Validation cohort 
Valpredgps<-predict(rf,data.test)
ConMat <- confusionMatrix(Valpredgps,pred.var)
ConMat
F1 <- ConMat[["byClass"]][["F1"]]
F1
Table <- as.vector(ConMat[["table"]])
Tablem <-matrix(c(Table),nrow=2)
MatCC <- mcc(confusionM=Tablem)
MatCC

#Validation AUC
Valvotes_Comb2=predict(rf,newdata = data.test,type = "prob")
write.csv((data.frame(Valvotes_Comb2)),file="Valvotes_Comb2F.csv")
Valpred=prediction(Valvotes_Comb2[,2],pred.var)
#Area under curve
Valauc=performance(Valpred, "auc")
Valauc <- Valauc@y.values[[1]]
Valauc 
#True Positive and Negative Rate
Valtpfp=performance(Valpred, "tpr","fpr")

#Results
write.csv((data.frame(OOBauc,Valauc,F1, MatCC)),file="Perf_Comb2F.csv")

OFPR <- unlist(OOBtpfp@x.values)
OTPR <- unlist(OOBtpfp@y.values)
OROC <- data.frame(OFPR,OTPR)
write.csv(OROC,file="OROC_Comb2F.csv")

VFPR <- unlist(Valtpfp@x.values)
VTPR <- unlist(Valtpfp@y.values)
VROC <- data.frame(VFPR,VTPR)
write.csv(VROC,file="VROC_Comb2F.csv")

summary(data.train$pred.var)
summary(pred.var)

#Boruta & RF - Comb3 (Parts 1,3,4): train, Part 2: test----
set.seed (123)
pred.var <- Comb3$Group
data.train<-cbind(Top500_Comb3,pred.var)

#Boruta
bor<-Boruta(pred.var~.,data=data.train,maxRuns=1000)
bor
plot(bor)
Boruta_Comb3_imphist <- bor[["ImpHistory"]]
write.csv(bor[["ImpHistory"]],file="Boruta_Comb3_imphistoryfullF.csv")
write.csv(bor[["finalDecision"]],file="Boruta_Comb3F.csv")
imp.probes_Comb3<-names(bor$finalDecision[which(bor$finalDecision=="Confirmed")])

data.train<-data.frame(cbind(Top500_Comb3[match(imp.probes_Comb3,colnames(Top500_Comb3))],pred.var))

#Running RF
rf<-randomForest(pred.var~.,data=data.train,ntree=1000)
rf

#OOB AUC
OOBvotes_Comb3=predict(rf,type = "prob")
write.csv((data.frame(OOBvotes_Comb3)),file="OOBvotes_Comb3F.csv")
OOBpred=prediction(OOBvotes_Comb3[,2],data.train$pred.var)
# Area under curve
OOBauc=performance(OOBpred, "auc")
OOBauc <- OOBauc@y.values[[1]]
OOBauc
# True Positive and Negative Rate
OOBtpfp=performance(OOBpred, "tpr","fpr")

#Preparation of test data (unsupervised)
pred.var<-as.factor(Part2$Group) #predictive variable
data.test<- Top500_Part2[,match(imp.probes_Comb3,colnames(Top500_Part2))]

#Validation cohort 
Valpredgps<-predict(rf,data.test)
ConMat <- confusionMatrix(Valpredgps,pred.var)
ConMat
F1 <- ConMat[["byClass"]][["F1"]]
F1
Table <- as.vector(ConMat[["table"]])
Tablem <-matrix(c(Table),nrow=2)
MatCC <- mcc(confusionM=Tablem)
MatCC

#Validation AUC
Valvotes_Comb3=predict(rf,newdata = data.test,type = "prob")
write.csv((data.frame(Valvotes_Comb3)),file="Valvotes_Comb3F.csv")
Valpred=prediction(Valvotes_Comb3[,2],pred.var)
#Area under curve
Valauc=performance(Valpred, "auc")
Valauc <- Valauc@y.values[[1]]
Valauc 
#True Positive and Negative Rate
Valtpfp=performance(Valpred, "tpr","fpr")

#Results
write.csv((data.frame(OOBauc,Valauc,F1, MatCC)),file="Perf_Comb3F.csv")

OFPR <- unlist(OOBtpfp@x.values)
OTPR <- unlist(OOBtpfp@y.values)
OROC <- data.frame(OFPR,OTPR)
write.csv(OROC,file="OROC_Comb3F.csv")

VFPR <- unlist(Valtpfp@x.values)
VTPR <- unlist(Valtpfp@y.values)
VROC <- data.frame(VFPR,VTPR)
write.csv(VROC,file="VROC_Comb3F.csv")

#Boruta & RF - Comb4 (Parts 2,3,4): train, Part 1: test----
set.seed (123)
pred.var <- Comb4$Group
data.train<-cbind(Top500_Comb4,pred.var)

#Boruta
bor<-Boruta(pred.var~.,data=data.train,maxRuns=1000)
bor
plot(bor)
Boruta_Comb4_imphist <- bor[["ImpHistory"]]
write.csv(bor[["ImpHistory"]],file="Boruta_Comb4_imphistoryfullF.csv")
write.csv(bor[["finalDecision"]],file="Boruta_Comb4F.csv")
imp.probes_Comb4<-names(bor$finalDecision[which(bor$finalDecision=="Confirmed")])

data.train<-data.frame(cbind(Top500_Comb4[match(imp.probes_Comb4,colnames(Top500_Comb4))],pred.var))

#Running RF
rf<-randomForest(pred.var~.,data=data.train,ntree=1000)
rf

#OOB AUC
OOBvotes_Comb4=predict(rf,type = "prob")
write.csv((data.frame(OOBvotes_Comb4)),file="OOBvotes_Comb4F.csv")
OOBpred=prediction(OOBvotes_Comb4[,2],data.train$pred.var)
# Area under curve
OOBauc=performance(OOBpred, "auc")
OOBauc <- OOBauc@y.values[[1]]
OOBauc
# True Positive and Negative Rate
OOBtpfp=performance(OOBpred, "tpr","fpr")

#Preparation of test data (unsupervised)
pred.var<-as.factor(Part1$Group) #predictive variable
data.test<- Top500_Part1[,match(imp.probes_Comb4,colnames(Top500_Part1))]

#Validation cohort 
Valpredgps<-predict(rf,data.test)
ConMat <- confusionMatrix(Valpredgps,pred.var)
ConMat
F1 <- ConMat[["byClass"]][["F1"]]
F1
Table <- as.vector(ConMat[["table"]])
Tablem <-matrix(c(Table),nrow=2)
MatCC <- mcc(confusionM=Tablem)
MatCC

#Validation AUC
Valvotes_Comb4=predict(rf,newdata = data.test,type = "prob")
write.csv((data.frame(Valvotes_Comb4)),file="Valvotes_Comb4F.csv")
Valpred=prediction(Valvotes_Comb4[,2],pred.var)
#Area under curve
Valauc=performance(Valpred, "auc")
Valauc <- Valauc@y.values[[1]]
Valauc 
#True Positive and Negative Rate
Valtpfp=performance(Valpred, "tpr","fpr")

#Results
write.csv((data.frame(OOBauc,Valauc,F1, MatCC)),file="Perf_Comb4F.csv")

OFPR <- unlist(OOBtpfp@x.values)
OTPR <- unlist(OOBtpfp@y.values)
OROC <- data.frame(OFPR,OTPR)
write.csv(OROC,file="OROC_Comb4F.csv")

VFPR <- unlist(Valtpfp@x.values)
VTPR <- unlist(Valtpfp@y.values)
VROC <- data.frame(VFPR,VTPR)
write.csv(VROC,file="VROC_Comb4F.csv")

#DMP info for predictive probes----
DMP_pred <- DMPtop500F[match(imp.probes_all,DMPtop500F$X),]
write.csv(DMP_pred, file="DMP_PredF.csv")

imp.probes_incdups <- c(imp.probes_Comb1,imp.probes_Comb2, imp.probes_Comb3, imp.probes_Comb4)
DMP_pred_incdups <- DMPtop500F[match(imp.probes_incdups,DMPtop500F$X),]
write.csv(DMP_pred_incdups, file="DMP_Pred_incdupsF.csv")