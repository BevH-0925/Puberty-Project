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
B2normc_top500 <- read.csv("~/Documents/Code/B2normc_top500.csv")
B2QMGp <- read.csv("~/Documents/Code/B2QMGp.csv")
load("~/Documents/Code/B2norm.Robj")
  
#Add/update rownames
rownames(B2QMGp)<-paste(B2QMGp[,1])
B2QMGp <- B2QMGp[,-1]

rownames(B2normc_top500)<-paste(B2normc_top500[,1])
B2normc_top500 <- B2normc_top500[,-1]

#Create a vector of DMP names (from top 500 Cl_B2)
DMP_top500_Cl_B2 <- colnames(B2normc_top500)

#SMOTE whole dataset
set.seed(123)
pred.var<-as.factor(B2QMGp$ClQGp) #predictive variable
data<-SMOTE(B2norm,pred.var,K=5)
data<-as.data.frame(data$data)
class<-as.factor(data$class) #predictive variable
data_1st <- data[data$class=="1st Quartile",]
data_2ndto4th <- data[data$class=="2nd-4th Quartiles",]
class1st <- as.factor(data_1st$class)
class2ndto4th <- as.factor(data_2ndto4th$class)
data_1st<- data_1st[,match(DMP_top500_Cl_B2,colnames(data_1st))]
data_2ndto4th<- data_2ndto4th[,match(DMP_top500_Cl_B2,colnames(data_2ndto4th))]
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
plot(bor, xlab="Top 500 probes for time to B2 (Q1 vs Q2-4)", ylab=NA, xaxt = "n", cex.lab=3, cex.axis=3)
mtext("Boruta Probe Importance Score", side=2, line=4, cex=3)
write.csv(bor[["finalDecision"]],file="Boruta_B2_cl_DMP.csv")

#Selecting genes confirmed to be important
imp.genes<-names(bor$finalDecision[which(bor$finalDecision=="Confirmed")])
write.csv(bor[["ImpHistory"]],file="ImpHistory_B2_cl_DMP.csv")
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
write.csv((data.frame(OOBauc,Valauc,F1,MatCC)),file="B2_cl_DMP.csv")

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
ggsave("ROC_B2_cl_DMP.png", dpi=320)

summary(data.train$pred.var)
summary(pred.var)

