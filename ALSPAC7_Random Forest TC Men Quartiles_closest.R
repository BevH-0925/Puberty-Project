library(randomForest)
library(Boruta)
library(pROC)
library(ROCR)
library(smotefamily)
library(caret)
library(mltools)
library(ggplot2)

setwd("~/Documents/Code/ALSPAC_Random Forest Analysis")

#Load data----
Expdata_top500Clmen <- read.csv("~/Documents/Code/Expdata_Cl_or.csv")
MenQgroups <- read.csv("~/Documents/Code/Puberty.DQG_or.csv")

rownames(Expdata_top500Clmen)<-paste(Expdata_top500Clmen[,1])
Expdata_top500Clmen <-Expdata_top500Clmen[,-1]

rownames(MenQgroups)<-paste(MenQgroups[,2])
MenQgroups <-MenQgroups[,-c(1:2)]

#Transcriptome data
load("~/Documents/Code/EXP.Rdata")

#Transpose transcriptomic data so individuals are in rows, genes in columns
Expdata_t <- (t(data))
rm(data)

#Remove row 642 and 728 due to -9; Patient IDs #-----A and ----A
Expdata_t <- Expdata_t[-c(642,728),]

#Subset of transcriptomic data to include only those participants with a time to menarche
Expdata_Men <- Expdata_t[match(rownames(MenQgroups),rownames(Expdata_t)),]
Expdata_Men <- as.data.frame(Expdata_Men)

#Create a vector of probe names (from top 500 Cl_Men)
Probes_top500_Cl_Men <- colnames(Expdata_top500Clmen)

#SMOTE whole dataset
set.seed(123)
pred.var<-as.factor(MenQgroups$ClQGp) #predictive variable
data<-SMOTE(Expdata_Men,pred.var,K=5)
data<-as.data.frame(data$data)
class<-as.factor(data$class) #predictive variable
data_1st <- data[data$class=="1st Quartile",]
data_2ndto4th <- data[data$class=="2nd-4th Quartiles",]
class1st <- as.factor(data_1st$class)
class2ndto4th <- as.factor(data_2ndto4th$class)
data_1st<- data_1st[,match(Probes_top500_Cl_Men,colnames(data_1st))]
data_2ndto4th<- data_2ndto4th[,match(Probes_top500_Cl_Men,colnames(data_2ndto4th))]
data_1st<-data.frame(cbind(data_1st,class1st))
data_2ndto4th<-data.frame(cbind(data_2ndto4th,class2ndto4th))

#1st & 2nd-4th split 70/30 into test & train
selected.samples1st<-sample(c(1:nrow(data_1st)),size = nrow(data_1st)*0.7)
train.data_1st<-data_1st[selected.samples1st,]
test.data_1st<-data_1st[-selected.samples1st,]

selected.samples2ndto4th<-sample(c(1:nrow(data_2ndto4th)),size = nrow(data_2ndto4th)*0.7)
train.data_2ndto4th<-data_2ndto4th[selected.samples2ndto4th,]
test.data_2ndto4th<-data_2ndto4th[-selected.samples2ndto4th,]
colnames(train.data_1st)[colnames(train.data_1st)=="class1st"] <- "class"
colnames(train.data_2ndto4th)[colnames(train.data_2ndto4th)=="class2ndto4th"] <- "class"
colnames(test.data_1st)[colnames(test.data_1st)=="class1st"] <- "class"
colnames(test.data_2ndto4th)[colnames(test.data_2ndto4th)=="class2ndto4th"] <- "class"
data.train<-rbind(train.data_1st,train.data_2ndto4th)
data.test<-rbind(test.data_1st,test.data_2ndto4th)

#Boruta
bor<-Boruta(class~.,data=data.train)
bor
par(mar = c(5, 7, 4, 2)) #Extra code to add a margin so y-axis label doesn't get cut-off
plot(bor, xlab="Top 500 probes for time to menarche (Q1 vs Q2-4)", ylab=NA, xaxt = "n", cex.lab=3, cex.axis=3)
mtext("Boruta Probe Importance Score", side=2, line=4, cex=3)
write.csv(bor[["finalDecision"]],file="Boruta_Men_cl_T.csv")

#Selecting genes confirmed to be important
imp.genes<-names(bor$finalDecision[which(bor$finalDecision=="Confirmed")])
write.csv(bor[["ImpHistory"]],file="ImpHistory_Men_cl_T.csv")
pred.var <- data.train$class
data.train<-data.frame(cbind(data.train[match(imp.genes,colnames(data.train))],pred.var))

#Running RF
rf<-randomForest(pred.var~.,data=data.train,ntree=1000)
rf

#OOB AUC
OOBvotes=predict(rf,type = "prob")
OOBpred=prediction(OOBvotes[,2],data.train$pred.var)
#Area under curve
OOBauc=performance(OOBpred, "auc")
OOBauc<- OOBauc@y.values[[1]]
OOBauc
#True Positive and Negative Rate
OOBtpfp=performance(OOBpred, "tpr","fpr")

#Preparation of test data (unsupervised)
pred.var<-as.factor(data.test$class) #predictive variable
data.test<- data.test[,match(imp.genes,colnames(data.test))]

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
Valpred=prediction(Valvotes[,2],pred.var)
#Area under curve
Valauc=performance(Valpred, "auc")
Valauc <- Valauc@y.values[[1]]
Valauc 
#True Positive and Negative Rate
Valtpfp=performance(Valpred, "tpr","fpr")

#Results
write.csv((data.frame(OOBauc,Valauc,F1,MatCC)),file="Men_cl_T.csv")

OFPR <- unlist(OOBtpfp@x.values)
OTPR <- unlist(OOBtpfp@y.values)
OROC<- data.frame(OFPR,OTPR)

VFPR <- unlist(Valtpfp@x.values)
VTPR <- unlist(Valtpfp@y.values)
VROC <- data.frame(VFPR,VTPR)

ggplot()+
  labs(x="False Positive Rate", y="True Positive Rate")+
  geom_line(data=OROC, aes(x=OFPR,y=OTPR,colour="Out of Bag"),show.legend=TRUE,lwd=1)+
  theme_classic(base_size = 18)+
  geom_line(data=VROC, aes(x=VFPR,y=VTPR,colour="Test Data"), show.legend=TRUE,lwd=1)+
  geom_abline(slope=1,intercept=0,lwd=1,lty=2,col="gray")+
  scale_color_discrete(name="Data Set", labels=c("Out of Bag (Training Data)", "Test Data Set"))
ggsave("ROC_Men_cl_T.png", dpi=320)

summary(data.train$pred.var)
summary(pred.var)

