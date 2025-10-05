library(randomForest)
library(Boruta)
library(pROC)
library(ROCR)
library(smotefamily)
library(caret)
library(mltools)
library(ggplot2)

setwd("~/Documents/RNASeqBH1-BH32/Random Forest")

#Prepare data----
Top500 <- read.csv("~/Documents/RNASeqBH1-BH32/EdgeR/Data/BH1-32T_top500_lcpm_classic.csv",row.names=1)

clindata <- read.csv("~/Documents/RNASeqBH1-BH32/Clinical Data Analysis/Data/BH1-32.csv", row.names=1)
clindata$Group<- as.factor(clindata$Group)

#Transpose
Top500 <- t(Top500)

#Divide into pubertal and pre-pubertal
clindata_prepub <- clindata[clindata$Group== "Prepubertal",]
clindata_pub <- clindata[clindata$Group== "Pubertal",]

#Randomly Split Dataset
set.seed(123) 
n <- nrow(clindata_prepub)  # Number of rows in the dataframe
split_ids <- sample(rep(1:4, length.out = n))  # Assign random numbers 1 to 4
clindata_prepub_parts <- split(clindata_prepub, split_ids)  # Split dataframe into 4 parts

n <- nrow(clindata_pub)  # Number of rows in the dataframe
split_ids <- sample(rep(1:4, length.out = n))  # Assign random numbers 1 to 4
clindata_pub_parts <- split(clindata_pub, split_ids)  # Split dataframe into 4 parts

Part1 <- rbind(clindata_prepub_parts[["1"]],clindata_pub_parts[["1"]])
Part2 <- rbind(clindata_prepub_parts[["2"]],clindata_pub_parts[["2"]])
Part3 <- rbind(clindata_prepub_parts[["3"]],clindata_pub_parts[["3"]])
Part4 <- rbind(clindata_prepub_parts[["4"]],clindata_pub_parts[["4"]])

Comb1 <- rbind(Part1,Part2,Part3)
Comb2 <- rbind(Part1,Part2,Part4)
Comb3 <- rbind(Part1,Part3,Part4)
Comb4 <- rbind(Part2,Part3,Part4)

Top500_Part1 <- as.data.frame(Top500[match(row.names(Part1),rownames(Top500)),])
Top500_Part2 <- as.data.frame(Top500[match(row.names(Part2),rownames(Top500)),])
Top500_Part3 <- as.data.frame(Top500[match(row.names(Part3),rownames(Top500)),])
Top500_Part4 <- as.data.frame(Top500[match(row.names(Part4),rownames(Top500)),])
Top500_Comb1 <- as.data.frame(Top500[match(row.names(Comb1),rownames(Top500)),])
Top500_Comb2 <- as.data.frame(Top500[match(row.names(Comb2),rownames(Top500)),])
Top500_Comb3 <- as.data.frame(Top500[match(row.names(Comb3),rownames(Top500)),])
Top500_Comb4 <- as.data.frame(Top500[match(row.names(Comb4),rownames(Top500)),])

#Boruta & RF - Comb1 (Parts 1,2,3): train, Part 4: test----
set.seed (123)
pred.var <- Comb1$Group
data.train<-cbind(Top500_Comb1,pred.var)

#Boruta
bor<-Boruta(pred.var~.,data=data.train,maxRuns=1000)
bor
plot(bor)

write.csv(bor[["finalDecision"]],file="Boruta_Comb1.csv")
imp.genes<-names(bor$finalDecision[which(bor$finalDecision=="Confirmed")])

data.train<-data.frame(cbind(Top500_Comb1[match(imp.genes,colnames(Top500_Comb1))],pred.var))

#Running RF
rf<-randomForest(pred.var~.,data=data.train,ntree=1000)
rf

#OOB AUC
OOBvotes=predict(rf,type = "prob")
write.csv((data.frame(OOBvotes)),file="OOBvotes_Comb1.csv")

OOBpred=prediction(OOBvotes[,2],data.train$pred.var)
# Area under curve
OOBauc=performance(OOBpred, "auc")
OOBauc <- OOBauc@y.values[[1]]
OOBauc
# True Positive and Negative Rate
OOBtpfp=performance(OOBpred, "tpr","fpr")

#Preparation of test data (unsupervised)
pred.var<-as.factor(Part4$Group) #predictive variable
data.test<- Top500_Part4[,match(imp.genes,colnames(Top500_Part4))]

#Issue with format of some of the gene names: need . instead of -
colnames(data.test)[colnames(data.test)=='RP11-538P18.2'] <- 'RP11.538P18.2'
#etc.

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
Valvotes=predict(rf,newdata = data.test,type = "prob")
write.csv((data.frame(Valvotes)),file="Valvotes_Comb1.csv")
Valpred=prediction(Valvotes[,2],pred.var)
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

ggplot()+
  labs(title="ROC Curve: Random Forest (with Boruta) Using Top 500 DEGs", subtitle="MF, classic, Comb 1",x="False Positive Rate", y="True Positive Rate")+
  geom_line(data=OROC, aes(x=OFPR,y=OTPR,colour="Out of Bag"),show.legend=TRUE,lwd=1)+
  theme_classic(base_size = 18)+
  geom_line(data=VROC, aes(x=VFPR,y=VTPR,colour="Test Data"), show.legend=TRUE,lwd=1)+
  geom_abline(slope=1,intercept=0,lwd=1,lty=2,col="gray")+
  scale_color_discrete(name="Data Set", labels=c("Out of Bag (Training Data)", "Test Data Set"))

summary(data.train$pred.var)
summary(pred.var)

#Compare peak LH to votes for puberty
PeakLHvsPubVotes <- clindata[,c(20,23,28)]
PeakLHvsPubVotes <- PeakLHvsPubVotes[match(rownames(Valvotes),rownames(PeakLHvsPubVotes)),]
PeakLHvsPubVotes <- cbind(PeakLHvsPubVotes,Valvotes[,2])
names(PeakLHvsPubVotes)[4] <- "Puberty Votes (%)"
PeakLHvsPubVotes <- PeakLHvsPubVotes[-6,]
ggplot(PeakLHvsPubVotes, aes(x=PeakLH,y=`Puberty Votes (%)`))+
  geom_point(aes(colour=Group, shape=Group))+
  labs(title="Relationship of Peak LH to RF Votes",subtitle="MF, classic - validation, Comb1", x="Peak LH Concentration (IU/L)", y="Puberty Votes (%)")+
  scale_x_continuous(breaks = c(0,10,20,30,40,50,60))+
  geom_abline(slope=0,intercept=0.5,lwd=1,lty=2,col="gray")+
  theme_classic(base_size = 18)

PeakLHvsPubVotesOOB <- clindata[,c(20,23,28)]
PeakLHvsPubVotesOOB <- PeakLHvsPubVotesOOB[match(rownames(OOBvotes),rownames(PeakLHvsPubVotesOOB)),]
PeakLHvsPubVotesOOB <- cbind(PeakLHvsPubVotesOOB,OOBvotes[,2])
names(PeakLHvsPubVotesOOB)[4] <- "Puberty Votes (%)"
PeakLHvsPubVotesOOB <- PeakLHvsPubVotesOOB[-15,]
ggplot(PeakLHvsPubVotesOOB, aes(x=PeakLH,y=`Puberty Votes (%)`))+
  geom_point(aes(colour=Group, shape=Group))+
  labs(title="Relationship of Peak LH to RF Votes",subtitle="MF, classic - OOB, Comb1", x="Peak LH Concentration (IU/L)", y="Puberty Votes (%)")+
  scale_x_continuous(breaks = c(0,10,20,30,40,50,60))+
  geom_abline(slope=0,intercept=0.5,lwd=1,lty=2,col="gray")+
  theme_classic(base_size = 18)

#Boruta & RF - Comb2 (Parts 1,2,4): train, Part 3: test----
set.seed (123)
pred.var <- Comb2$Group
data.train<-cbind(Top500_Comb2,pred.var)

#Boruta
bor<-Boruta(pred.var~.,data=data.train,maxRuns=1000)
bor
plot(bor)

write.csv(bor[["finalDecision"]],file="Boruta_Comb2.csv")
imp.genes<-names(bor$finalDecision[which(bor$finalDecision=="Confirmed")])

data.train<-data.frame(cbind(Top500_Comb2[match(imp.genes,colnames(Top500_Comb2))],pred.var))

#Running RF
rf<-randomForest(pred.var~.,data=data.train,ntree=1000)
rf

#OOB AUC
OOBvotes=predict(rf,type = "prob")
write.csv((data.frame(OOBvotes)),file="OOBvotes_Comb2.csv")
OOBpred=prediction(OOBvotes[,2],data.train$pred.var)
# Area under curve
OOBauc=performance(OOBpred, "auc")
OOBauc <- OOBauc@y.values[[1]]
OOBauc
# True Positive and Negative Rate
OOBtpfp=performance(OOBpred, "tpr","fpr")

#Preparation of test data (unsupervised)
pred.var<-as.factor(Part3$Group) #predictive variable
data.test<- Top500_Part3[,match(imp.genes,colnames(Top500_Part3))]

#Issue with format of some of the gene names: need . instead of -
colnames(data.test)[colnames(data.test)=='RP11-538P18.2'] <- 'RP11.538P18.2'
#etc.

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
Valvotes=predict(rf,newdata = data.test,type = "prob")
write.csv((data.frame(Valvotes)),file="Valvotes_Comb2.csv")
Valpred=prediction(Valvotes[,2],pred.var)
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

ggplot()+
  labs(title="ROC Curve: Random Forest (with Boruta) Using Top 500 DEGs", subtitle="MF, classic, Comb 2",x="False Positive Rate", y="True Positive Rate")+
  geom_line(data=OROC, aes(x=OFPR,y=OTPR,colour="Out of Bag"),show.legend=TRUE,lwd=1)+
  theme_classic(base_size = 18)+
  geom_line(data=VROC, aes(x=VFPR,y=VTPR,colour="Test Data"), show.legend=TRUE,lwd=1)+
  geom_abline(slope=1,intercept=0,lwd=1,lty=2,col="gray")+
  scale_color_discrete(name="Data Set", labels=c("Out of Bag (Training Data)", "Test Data Set"))

summary(data.train$pred.var)
summary(pred.var)

#Compare peak LH to votes for puberty
PeakLHvsPubVotes <- clindata[,c(20,23,28)]
PeakLHvsPubVotes <- PeakLHvsPubVotes[match(rownames(Valvotes),rownames(PeakLHvsPubVotes)),]
PeakLHvsPubVotes <- cbind(PeakLHvsPubVotes,Valvotes[,2])
names(PeakLHvsPubVotes)[4] <- "Puberty Votes (%)"
PeakLHvsPubVotes <- PeakLHvsPubVotes[-6,]
ggplot(PeakLHvsPubVotes, aes(x=PeakLH,y=`Puberty Votes (%)`))+
  geom_point(aes(colour=Group, shape=Group))+
  labs(title="Relationship of Peak LH to RF Votes",subtitle="MF, classic - validation, Comb2", x="Peak LH Concentration (IU/L)", y="Puberty Votes (%)")+
  scale_x_continuous(breaks = c(0,10,20,30,40,50,60))+
  geom_abline(slope=0,intercept=0.5,lwd=1,lty=2,col="gray")+
  theme_classic(base_size = 18)

PeakLHvsPubVotesOOB <- clindata[,c(20,23,28)]
PeakLHvsPubVotesOOB <- PeakLHvsPubVotesOOB[match(rownames(OOBvotes),rownames(PeakLHvsPubVotesOOB)),]
PeakLHvsPubVotesOOB <- cbind(PeakLHvsPubVotesOOB,OOBvotes[,2])
names(PeakLHvsPubVotesOOB)[4] <- "Puberty Votes (%)"
PeakLHvsPubVotesOOB <- PeakLHvsPubVotesOOB[-15,]
ggplot(PeakLHvsPubVotesOOB, aes(x=PeakLH,y=`Puberty Votes (%)`))+
  geom_point(aes(colour=Group, shape=Group))+
  labs(title="Relationship of Peak LH to RF Votes",subtitle="MF, classic - OOB, Comb2", x="Peak LH Concentration (IU/L)", y="Puberty Votes (%)")+
  scale_x_continuous(breaks = c(0,10,20,30,40,50,60))+
  geom_abline(slope=0,intercept=0.5,lwd=1,lty=2,col="gray")+
  theme_classic(base_size = 18)

#Boruta & RF - Comb3 (Parts 1,3,4): train, Part 2: test----
set.seed (123)
pred.var <- Comb3$Group
data.train<-cbind(Top500_Comb3,pred.var)

#Boruta
bor<-Boruta(pred.var~.,data=data.train,maxRuns=1000)
bor
plot(bor)

write.csv(bor[["finalDecision"]],file="Boruta_Comb3.csv")
imp.genes<-names(bor$finalDecision[which(bor$finalDecision=="Confirmed")])

data.train<-data.frame(cbind(Top500_Comb3[match(imp.genes,colnames(Top500_Comb3))],pred.var))

#Running RF
rf<-randomForest(pred.var~.,data=data.train,ntree=1000)
rf

#OOB AUC
OOBvotes=predict(rf,type = "prob")
write.csv((data.frame(OOBvotes)),file="OOBvotes_Comb3.csv")
OOBpred=prediction(OOBvotes[,2],data.train$pred.var)
# Area under curve
OOBauc=performance(OOBpred, "auc")
OOBauc <- OOBauc@y.values[[1]]
OOBauc
# True Positive and Negative Rate
OOBtpfp=performance(OOBpred, "tpr","fpr")

#Preparation of test data (unsupervised)
pred.var<-as.factor(Part2$Group) #predictive variable
data.test<- Top500_Part2[,match(imp.genes,colnames(Top500_Part2))]

#Issue with format of some of the gene names: need . instead of -
colnames(data.test)[colnames(data.test)=='RP11-538P18.2'] <- 'RP11.538P18.2'
#etc.

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
Valvotes=predict(rf,newdata = data.test,type = "prob")
write.csv((data.frame(Valvotes)),file="Valvotes_Comb3.csv")
Valpred=prediction(Valvotes[,2],pred.var)
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

ggplot()+
  labs(title="ROC Curve: Random Forest (with Boruta) Using Top 500 DEGs", subtitle="MF, classic, Comb 3",x="False Positive Rate", y="True Positive Rate")+
  geom_line(data=OROC, aes(x=OFPR,y=OTPR,colour="Out of Bag"),show.legend=TRUE,lwd=1)+
  theme_classic(base_size = 18)+
  geom_line(data=VROC, aes(x=VFPR,y=VTPR,colour="Test Data"), show.legend=TRUE,lwd=1)+
  geom_abline(slope=1,intercept=0,lwd=1,lty=2,col="gray")+
  scale_color_discrete(name="Data Set", labels=c("Out of Bag (Training Data)", "Test Data Set"))

summary(data.train$pred.var)
summary(pred.var)

#Compare peak LH to votes for puberty
PeakLHvsPubVotes <- clindata[,c(20,23,28)]
PeakLHvsPubVotes <- PeakLHvsPubVotes[match(rownames(Valvotes),rownames(PeakLHvsPubVotes)),]
PeakLHvsPubVotes <- cbind(PeakLHvsPubVotes,Valvotes[,2])
names(PeakLHvsPubVotes)[4] <- "Puberty Votes (%)"
PeakLHvsPubVotes <- PeakLHvsPubVotes[-6,]
ggplot(PeakLHvsPubVotes, aes(x=PeakLH,y=`Puberty Votes (%)`))+
  geom_point(aes(colour=Group, shape=Group))+
  labs(title="Relationship of Peak LH to RF Votes",subtitle="MF, classic - validation, Comb3", x="Peak LH Concentration (IU/L)", y="Puberty Votes (%)")+
  scale_x_continuous(breaks = c(0,10,20,30,40,50,60))+
  geom_abline(slope=0,intercept=0.5,lwd=1,lty=2,col="gray")+
  theme_classic(base_size = 18)

PeakLHvsPubVotesOOB <- clindata[,c(20,23,28)]
PeakLHvsPubVotesOOB <- PeakLHvsPubVotesOOB[match(rownames(OOBvotes),rownames(PeakLHvsPubVotesOOB)),]
PeakLHvsPubVotesOOB <- cbind(PeakLHvsPubVotesOOB,OOBvotes[,2])
names(PeakLHvsPubVotesOOB)[4] <- "Puberty Votes (%)"
PeakLHvsPubVotesOOB <- PeakLHvsPubVotesOOB[-15,]
ggplot(PeakLHvsPubVotesOOB, aes(x=PeakLH,y=`Puberty Votes (%)`))+
  geom_point(aes(colour=Group, shape=Group))+
  labs(title="Relationship of Peak LH to RF Votes",subtitle="MF, classic - OOB, Comb3", x="Peak LH Concentration (IU/L)", y="Puberty Votes (%)")+
  scale_x_continuous(breaks = c(0,10,20,30,40,50,60))+
  geom_abline(slope=0,intercept=0.5,lwd=1,lty=2,col="gray")+
  theme_classic(base_size = 18)

#Boruta & RF - Comb4 (Parts 2,3,4): train, Part 1: test----
set.seed (123)
pred.var <- Comb4$Group
data.train<-cbind(Top500_Comb4,pred.var)

#Boruta
bor<-Boruta(pred.var~.,data=data.train,maxRuns=1000)
bor
plot(bor)

write.csv(bor[["finalDecision"]],file="Boruta_Comb4.csv")
imp.genes<-names(bor$finalDecision[which(bor$finalDecision=="Confirmed")])

data.train<-data.frame(cbind(Top500_Comb4[match(imp.genes,colnames(Top500_Comb4))],pred.var))

#Running RF
rf<-randomForest(pred.var~.,data=data.train,ntree=1000)
rf

#OOB AUC
OOBvotes=predict(rf,type = "prob")
write.csv((data.frame(OOBvotes)),file="OOBvotes_Comb4.csv")
OOBpred=prediction(OOBvotes[,2],data.train$pred.var)
# Area under curve
OOBauc=performance(OOBpred, "auc")
OOBauc <- OOBauc@y.values[[1]]
OOBauc
# True Positive and Negative Rate
OOBtpfp=performance(OOBpred, "tpr","fpr")

#Preparation of test data (unsupervised)
pred.var<-as.factor(Part1$Group) #predictive variable
data.test<- Top500_Part1[,match(imp.genes,colnames(Top500_Part1))]

#Issue with format of some of the gene names: need . instead of -
colnames(data.test)[colnames(data.test)=='RP11-538P18.2'] <- 'RP11.538P18.2'
#etc.

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
Valvotes=predict(rf,newdata = data.test,type = "prob")
write.csv((data.frame(Valvotes)),file="Valvotes_Comb4.csv")
Valpred=prediction(Valvotes[,2],pred.var)
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

ggplot()+
  labs(title="ROC Curve: Random Forest (with Boruta) Using Top 500 DEGs", subtitle="MF, classic, Comb 4",x="False Positive Rate", y="True Positive Rate")+
  geom_line(data=OROC, aes(x=OFPR,y=OTPR,colour="Out of Bag"),show.legend=TRUE,lwd=1)+
  theme_classic(base_size = 18)+
  geom_line(data=VROC, aes(x=VFPR,y=VTPR,colour="Test Data"), show.legend=TRUE,lwd=1)+
  geom_abline(slope=1,intercept=0,lwd=1,lty=2,col="gray")+
  scale_color_discrete(name="Data Set", labels=c("Out of Bag (Training Data)", "Test Data Set"))

summary(data.train$pred.var)
summary(pred.var)

#Compare peak LH to votes for puberty
PeakLHvsPubVotes <- clindata[,c(20,23,28)]
PeakLHvsPubVotes <- PeakLHvsPubVotes[match(rownames(Valvotes),rownames(PeakLHvsPubVotes)),]
PeakLHvsPubVotes <- cbind(PeakLHvsPubVotes,Valvotes[,2])
names(PeakLHvsPubVotes)[4] <- "Puberty Votes (%)"
PeakLHvsPubVotes <- PeakLHvsPubVotes[-6,]
ggplot(PeakLHvsPubVotes, aes(x=PeakLH,y=`Puberty Votes (%)`))+
  geom_point(aes(colour=Group, shape=Group))+
  labs(title="Relationship of Peak LH to RF Votes",subtitle="MF, classic - validation, Comb4", x="Peak LH Concentration (IU/L)", y="Puberty Votes (%)")+
  scale_x_continuous(breaks = c(0,10,20,30,40,50,60))+
  geom_abline(slope=0,intercept=0.5,lwd=1,lty=2,col="gray")+
  theme_classic(base_size = 18)

PeakLHvsPubVotesOOB <- clindata[,c(20,23,28)]
PeakLHvsPubVotesOOB <- PeakLHvsPubVotesOOB[match(rownames(OOBvotes),rownames(PeakLHvsPubVotesOOB)),]
PeakLHvsPubVotesOOB <- cbind(PeakLHvsPubVotesOOB,OOBvotes[,2])
names(PeakLHvsPubVotesOOB)[4] <- "Puberty Votes (%)"
PeakLHvsPubVotesOOB <- PeakLHvsPubVotesOOB[-15,]
ggplot(PeakLHvsPubVotesOOB, aes(x=PeakLH,y=`Puberty Votes (%)`))+
  geom_point(aes(colour=Group, shape=Group))+
  labs(title="Relationship of Peak LH to RF Votes",subtitle="MF, classic - OOB, Comb4", x="Peak LH Concentration (IU/L)", y="Puberty Votes (%)")+
  scale_x_continuous(breaks = c(0,10,20,30,40,50,60))+
  geom_abline(slope=0,intercept=0.5,lwd=1,lty=2,col="gray")+
  theme_classic(base_size = 18)