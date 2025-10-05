setwd("~/Documents/Code/RMCH_cross_validation")

library(caret)
library(randomForest)
library(kernelshap)
library(shapviz)
library(pROC)
library(tidyverse)
library(ggrepel)

set.seed(123)

#RMCH Transcriptome MF, top 500 GLM----

#Prepare data
Top500 <- read.csv("~/Documents/Code/RNASeqBH1-BH32/EdgeR/Data/BH1-32T_top500_lcpm_glmsea.csv",row.names=1)
clindata <- read.csv("~/Documents/Code/RNASeqBH1-BH32/Clinical Data Analysis/Data/BH1-32.csv", row.names=1)

#Transpose
Top500 <- t(Top500)
Top500 <- as.data.frame(Top500)

#Ensure clinical data and transcriptome data in same order
clindata <- clindata[rownames(Top500), ]

Class <- clindata$Group
Class <- as.factor(Class)

dataT<-cbind(Top500,Class)

#Define stratified k-fold cross-validation
ctrl <- trainControl(method="cv", number=5,classProbs = TRUE,summaryFunction = twoClassSummary, savePredictions = "final")

#Train Random Forest model with stratified CV
set.seed(123)
rf_model <- train(Class~.,data=dataT, method="rf", trControl=ctrl, metric="ROC", tuneLength=5)
print(rf_model)
plot(rf_model)

#Calculate OOB AUC
final_rf <- rf_model$finalModel

# OOB predictions (votes or probabilities) and true classes
oob_pred <- final_rf$predicted  # factor predictions
oob_prob <- final_rf$votes[,2]  # probability for positive class
oob_true <- rf_model$trainingData$.outcome  # true classes in training set

# Calculate OOB AUC
roc_obj <- roc(oob_true, oob_prob)
oob_auc <- auc(roc_obj)
print(oob_auc)

#RMCH Methylome MF top 500 DMPs----
clindataMeth <- read.csv("~/Documents/Code/RMCH_MethArrayBatch1/RMCHMethChAMP/Data/Meth1to46_prelim.csv")
DMPtop500list <- read.csv("~/Documents/Code/RMCH_MethArrayBatch1/RMCHMethChAMP/DMPtop500.csv")
load("~/Documents/Code/RMCH_MethArrayBatch1/RMCHMethChAMP/MyNormt.Robj")

matchingCpGstop500 <- intersect(names(MyNormt), DMPtop500list$Probe)
RMCH_M_top500_n46<- MyNormt[,matchingCpGstop500]

#Ensure clinical data and methylome data in same order
rownames(clindataMeth) <- clindataMeth$ParticipantID
clindataMeth <- clindataMeth[rownames(RMCH_M_top500_n46), ]

Class <- clindataMeth$Group
Class <- as.factor(Class)
RMCH_M_top500_n46 <- as.data.frame(RMCH_M_top500_n46)
dataM<-cbind(RMCH_M_top500_n46,Class)

#Train Random Forest model with stratified CV
set.seed(123)
rf_modelM <- train(Class~.,data=dataM, method="rf", trControl=ctrl, metric="ROC", tuneLength=5)

print(rf_modelM)
plot(rf_modelM)

#write.csv(dataM,"RMCH_M_top500_n46.csv")

#Calculate OOB AUC
final_rf <- rf_modelM$finalModel

# OOB predictions (votes or probabilities) and true classes
oob_pred <- final_rf$predicted  # factor predictions
oob_prob <- final_rf$votes[,2]  # probability for positive class
oob_true <- rf_modelM$trainingData$.outcome  # true classes in training set

# Calculate OOB AUC
roc_obj <- roc(oob_true, oob_prob)
oob_auc <- auc(roc_obj)
print(oob_auc)

#Repeat with ALSPAC predictive genes - transcriptome----
RMCH_T_matchingALSPACPredGenes <- read.csv("~/Documents/Code/ALSPAC and RMCH/ALSPAC model with genes common to RNASeq/RMCH_T_matchingALSPACPredGenes.csv")
rownames(RMCH_T_matchingALSPACPredGenes) <- RMCH_T_matchingALSPACPredGenes$X
RMCH_T_matchingALSPACPredGenes <- RMCH_T_matchingALSPACPredGenes[,-1]

#Ensure clinical data and transcriptome data in same order
clindata <- clindata[rownames(RMCH_T_matchingALSPACPredGenes), ]

Class <- clindata$Group
Class <- as.factor(Class)

dataTP<-cbind(RMCH_T_matchingALSPACPredGenes,Class)

#Define stratified k-fold cross-validation
ctrl <- trainControl(method="cv", number=5,classProbs = TRUE,summaryFunction = twoClassSummary, savePredictions = "final")

#Train Random Forest model with stratified CV
set.seed(123)
rf_modelP <- train(Class~.,data=dataTP, method="rf", trControl=ctrl, metric="ROC", tuneLength=5)

print(rf_modelP)
plot(rf_modelP)

resampled <- rf_modelP[["resampledCM"]]
pred <- rf_modelP[["pred"]]

#Calculate OOB AUC
final_rf <- rf_modelP$finalModel

# OOB predictions (votes or probabilities) and true classes
oob_pred <- final_rf$predicted  # factor predictions
oob_prob <- final_rf$votes[,2]  # probability for positive class
oob_true <- rf_modelP$trainingData$.outcome  # true classes in training set

# Calculate OOB AUC
roc_obj <- roc(oob_true, oob_prob)
oob_auc <- auc(roc_obj)
print(oob_auc)

#Repeat with ALSPAC predictive CpGs - methylome----
clindataMeth <- read.csv("~/Documents/Code/RMCH_MethArrayBatch1/RMCHMethChAMP/Data/Meth1to46_prelim.csv")
RMCHMethmatchingALSPACPredCpGs <- read.csv("~/Documents/Code/ALSPAC and RMCH/Methylome ALSPAC and RMCH/RMCHMethmatchingALSPACPredCpGs.csv")
rownames(clindataMeth) <- clindataMeth$ParticipantID
rownames(RMCHMethmatchingALSPACPredCpGs) <- RMCHMethmatchingALSPACPredCpGs$X
RMCHMethmatchingALSPACPredCpGs <- RMCHMethmatchingALSPACPredCpGs[,-1]

#Ensure clinical data and methylome data in same order
clindataMeth <- clindataMeth[rownames(RMCHMethmatchingALSPACPredCpGs), ]

Class <- clindataMeth$Group
Class <- as.factor(Class)
dataMP<-cbind(RMCHMethmatchingALSPACPredCpGs,Class)

#Train Random Forest model with stratified CV
set.seed(123)
rf_modelMP <- train(Class~.,data=dataMP, method="rf", trControl=ctrl, metric="ROC", tuneLength=5)

print(rf_modelMP)
plot(rf_modelMP)

#Calculate OOB AUC
final_rf <- rf_modelMP$finalModel

# OOB predictions (votes or probabilities) and true classes
oob_pred <- final_rf$predicted  # factor predictions
oob_prob <- final_rf$votes[,2]  # probability for positive class
oob_true <- rf_modelMP$trainingData$.outcome  # true classes in training set

# Calculate OOB AUC
roc_obj <- roc(oob_true, oob_prob)
oob_auc <- auc(roc_obj)
print(oob_auc)

resampled <- rf_modelMP[["resampledCM"]]
write.csv(resampled, "FoldperformanceRMCHMALSPACgenes.csv")
pred <- rf_modelMP[["pred"]]
write.csv(pred, "PredRMCHMALSPACgenes.csv")

#Shapley Values----
#RMCH Transcriptome, ALSPAC predictive genes
set.seed(123)
s <- kernelshap(rf_modelP, dataTP[-157], bg_X=dataTP, type="prob")

sv <- shapviz(s)
png(filename = "SHAPPlot_RMCHT_ALSPACPredgenes_320dpi.png", width = 6 * 320, height = 4 * 320, res = 320)
sv_importance(sv, kind = "bee")
dev.off()

#RMCH Methylome, ALSPAC predictive genes
set.seed(123)
sM <- kernelshap(rf_modelMP, dataMP[-386], bg_X=dataMP, type="prob")
svM <- shapviz(sM)
png(filename = "SHAPPlot_RMCHM_ALSPACPredCpGs_320dpi.png", width = 6 * 320, height = 4 * 320, res = 320)
sv_importance(svM, kind = "bee")
dev.off()

#Create df of top 10 pred genes - RMCH T, ALPAC T Pred genes
Top10T <- dataTP[,c("DNMT3B", "TTLL11", "TMEM192", "WDR66", "COQ6", "PC", "COLQ", "IL21R", "DYRK1A", "PIH1D1", "Class")]
write.csv(Top10T, "Top10T.csv")

#Create df of top 10 pred CpGs - RMCH M, ALPAC M Pred genes
Top10M <- dataMP[,c("cg24387126", "cg27072813", "cg15626881", "cg22580353", "cg09877009", "cg08505473", "cg18984282", "cg15318627", "cg17421241", "cg25603742", "Class")]
write.csv(Top10M, "Top10M.csv")

#Logistic regression - Transcriptome----
# Fit logistic regression model for DNMT3B
model <- glm(Class ~ DNMT3B, data = Top10T, family = "binomial")

# Display model summary
summary(model)
oddsratio <- exp(coef(model))  # Converts log-odds coefficients to odds ratios
write.csv(oddsratio,"OR_dfDNMT3B.csv")

# Get the summary object
model_summary <- summary(model)

# Extract coefficients as a matrix
coeff_matrix <- model_summary$coefficients

# Convert to data frame
coeff_df <- as.data.frame(coeff_matrix)

# Optionally, add row names (term names) as a column
coeff_df$Term <- rownames(coeff_df)
rownames(coeff_df) <- NULL

# Reorder columns so Term is first
coeff_df <- coeff_df[, c(ncol(coeff_df), 1:(ncol(coeff_df)-1))]

# Print the data frame
print(coeff_df)

write.csv(coeff_df,"coeff_dfDNMT3B.csv")

predicted_probs <- predict(model, type = "response")
predicted_classes <- ifelse(predicted_probs > 0.5, 1, 0)

# Confusion matrix
table(Predicted = predicted_classes, Actual = Top10T$Class)

# Generate predicted probabilities from the model
predicted_probs <- predict(model, type = "response")

# # Create the ROC curve
roc_curve <- roc(Top10T$Class, predicted_probs)
auc_value <- auc(roc_curve)
auc_value

# Fit logistic regression model for DNMT3B, TTLL11
model <- glm(Class ~ DNMT3B + TTLL11, data = Top10T, family = "binomial")
summary(model)
oddsratio <- exp(coef(model)) 
write.csv(oddsratio,"OR_df2T.csv")

model_summary <- summary(model)
coeff_matrix <- model_summary$coefficients
coeff_df <- as.data.frame(coeff_matrix)
coeff_df$Term <- rownames(coeff_df)
rownames(coeff_df) <- NULL
coeff_df <- coeff_df[, c(ncol(coeff_df), 1:(ncol(coeff_df)-1))]
print(coeff_df)
write.csv(coeff_df,"coeff_df2T.csv")

predicted_probs <- predict(model, type = "response")
predicted_classes <- ifelse(predicted_probs > 0.5, 1, 0)

table(Predicted = predicted_classes, Actual = Top10T$Class)
predicted_probs <- predict(model, type = "response")
roc_curve <- roc(Top10T$Class, predicted_probs)
auc_value <- auc(roc_curve)
auc_value

# Fit logistic regression model for DNMT3B, TTLL11, TMEM192
model <- glm(Class ~ DNMT3B + TTLL11 + TMEM192, data = Top10T, family = "binomial")
summary(model)
oddsratio <- exp(coef(model)) 
write.csv(oddsratio,"OR_df3T.csv")

model_summary <- summary(model)
coeff_matrix <- model_summary$coefficients
coeff_df <- as.data.frame(coeff_matrix)
coeff_df$Term <- rownames(coeff_df)
rownames(coeff_df) <- NULL
coeff_df <- coeff_df[, c(ncol(coeff_df), 1:(ncol(coeff_df)-1))]
print(coeff_df)
write.csv(coeff_df,"coeff_df3T.csv")

predicted_probs <- predict(model, type = "response")
predicted_classes <- ifelse(predicted_probs > 0.5, 1, 0)

table(Predicted = predicted_classes, Actual = Top10T$Class)
predicted_probs <- predict(model, type = "response")
roc_curve <- roc(Top10T$Class, predicted_probs)
auc_value <- auc(roc_curve)
auc_value

# Fit logistic regression model for DNMT3B, TTLL11, TMEM192, WDR66
model <- glm(Class ~ DNMT3B + TTLL11 + TMEM192 + WDR66, data = Top10T, family = "binomial")
summary(model)
oddsratio <- exp(coef(model)) 
write.csv(oddsratio,"OR_df4T.csv") 
model_summary <- summary(model)
coeff_matrix <- model_summary$coefficients
coeff_df <- as.data.frame(coeff_matrix)
coeff_df$Term <- rownames(coeff_df)
rownames(coeff_df) <- NULL
coeff_df <- coeff_df[, c(ncol(coeff_df), 1:(ncol(coeff_df)-1))]
print(coeff_df)
write.csv(coeff_df,"coeff_df4T.csv")
predicted_probs <- predict(model, type = "response")
predicted_classes <- ifelse(predicted_probs > 0.5, 1, 0)
table(Predicted = predicted_classes, Actual = Top10T$Class)
predicted_probs <- predict(model, type = "response")
roc_curve <- roc(Top10T$Class, predicted_probs)
auc_value <- auc(roc_curve)
auc_value

# Fit logistic regression model for DNMT3B, TTLL11, TMEM192, WDR66, COQ6
model <- glm(Class ~ DNMT3B + TTLL11 + TMEM192 + WDR66 + COQ6, data = Top10T, family = "binomial")
summary(model)
oddsratio <- exp(coef(model)) 
write.csv(oddsratio,"OR_df5T.csv")
model_summary <- summary(model)
coeff_matrix <- model_summary$coefficients
coeff_df <- as.data.frame(coeff_matrix)
coeff_df$Term <- rownames(coeff_df)
rownames(coeff_df) <- NULL
coeff_df <- coeff_df[, c(ncol(coeff_df), 1:(ncol(coeff_df)-1))]
print(coeff_df)
write.csv(coeff_df,"coeff_df5T.csv")
predicted_probs <- predict(model, type = "response")
predicted_classes <- ifelse(predicted_probs > 0.5, 1, 0)
table(Predicted = predicted_classes, Actual = Top10T$Class)
predicted_probs <- predict(model, type = "response")
roc_curve <- roc(Top10T$Class, predicted_probs)
auc_value <- auc(roc_curve)
auc_value

# Fit logistic regression model for DNMT3B, TTLL11, TMEM192, WDR66, COQ6, PC
model <- glm(Class ~ DNMT3B + TTLL11 + TMEM192 + WDR66 + COQ6 + PC, data = Top10T, family = "binomial")
summary(model)
oddsratio <- exp(coef(model)) 
write.csv(oddsratio,"OR_df6T.csv")
model_summary <- summary(model)
coeff_matrix <- model_summary$coefficients
coeff_df <- as.data.frame(coeff_matrix)
coeff_df$Term <- rownames(coeff_df)
rownames(coeff_df) <- NULL
coeff_df <- coeff_df[, c(ncol(coeff_df), 1:(ncol(coeff_df)-1))]
print(coeff_df)
write.csv(coeff_df,"coeff_df6T.csv")
predicted_probs <- predict(model, type = "response")
predicted_classes <- ifelse(predicted_probs > 0.5, 1, 0)
table(Predicted = predicted_classes, Actual = Top10T$Class)
predicted_probs <- predict(model, type = "response")
roc_curve <- roc(Top10T$Class, predicted_probs)
auc_value <- auc(roc_curve)
auc_value

# Fit logistic regression model for DNMT3B, TTLL11, TMEM192, WDR66, COQ6, PC, COLQ
model <- glm(Class ~ DNMT3B + TTLL11 + TMEM192 + WDR66 + COQ6 + PC + COLQ, data = Top10T, family = "binomial")
summary(model)
oddsratio <- exp(coef(model)) 
write.csv(oddsratio,"OR_df7T.csv")
model_summary <- summary(model)
coeff_matrix <- model_summary$coefficients
coeff_df <- as.data.frame(coeff_matrix)
coeff_df$Term <- rownames(coeff_df)
rownames(coeff_df) <- NULL
coeff_df <- coeff_df[, c(ncol(coeff_df), 1:(ncol(coeff_df)-1))]
print(coeff_df)
write.csv(coeff_df,"coeff_df7T.csv")
predicted_probs <- predict(model, type = "response")
predicted_classes <- ifelse(predicted_probs > 0.5, 1, 0)
table(Predicted = predicted_classes, Actual = Top10T$Class)
predicted_probs <- predict(model, type = "response")
roc_curve <- roc(Top10T$Class, predicted_probs)
auc_value <- auc(roc_curve)
auc_value

# Fit logistic regression model for DNMT3B, TTLL11, TMEM192, WDR66, COQ6, PC, COLQ, IL21R
model <- glm(Class ~ DNMT3B + TTLL11 + TMEM192 + WDR66 + COQ6 + PC + COLQ + IL21R, data = Top10T, family = "binomial")
summary(model)
oddsratio <- exp(coef(model)) 
write.csv(oddsratio,"OR_df8T.csv")
model_summary <- summary(model)
coeff_matrix <- model_summary$coefficients
coeff_df <- as.data.frame(coeff_matrix)
coeff_df$Term <- rownames(coeff_df)
rownames(coeff_df) <- NULL
coeff_df <- coeff_df[, c(ncol(coeff_df), 1:(ncol(coeff_df)-1))]
print(coeff_df)
write.csv(coeff_df,"coeff_df8T.csv")
predicted_probs <- predict(model, type = "response")
predicted_classes <- ifelse(predicted_probs > 0.5, 1, 0)
table(Predicted = predicted_classes, Actual = Top10T$Class)
predicted_probs <- predict(model, type = "response")
roc_curve <- roc(Top10T$Class, predicted_probs)
auc_value <- auc(roc_curve)
auc_value


library(logistf)
fit <- logistf(Class ~ DNMT3B + TTLL11 + TMEM192 + WDR66 + COQ6 + PC + COLQ + IL21R, data = Top10T)
summary(fit)

oddsratio <- exp(coef(fit)) 
write.csv(oddsratio,"OR_df8T_logistF.csv")
model_summary <- summary(fit)
coeff_matrix <- model_summary$coefficients
coeff_df <- as.data.frame(coeff_matrix)
coeff_df$Term <- rownames(coeff_df)
rownames(coeff_df) <- NULL
coeff_df <- coeff_df[, c(ncol(coeff_df), 1:(ncol(coeff_df)-1))]
print(coeff_df)
write.csv(coeff_df,"coeff_df8TlogistF.csv")
pval <- model_summary$prob
write.csv(pval,"pval_df8TlogistF.csv")
# variance-covariance matrix of coefficients
var_cov_matrix <- fit$var

# standard errors are square root of diagonal of var covariance matrix
se <- sqrt(diag(var_cov_matrix))
write.csv(se,"se_df8TlogistF.csv")
CI <- exp(cbind(OR = coef(fit), confint(fit)))
write.csv(CI,"ORandCI_df8TlogistF.csv")

predicted_probs <- predict(fit, type = "response")
predicted_classes <- ifelse(predicted_probs > 0.5, 1, 0)
table(Predicted = predicted_classes, Actual = Top10T$Class)
predicted_probs <- predict(fit, type = "response")
roc_curve <- roc(Top10T$Class, predicted_probs)
auc_value <- auc(roc_curve)
auc_value

# Fit logistic regression model for DNMT3B, TTLL11, TMEM192, WDR66, COQ6, PC, COLQ, IL21R, DYRK1A
model <- glm(Class ~ DNMT3B + TTLL11 + TMEM192 + WDR66 + COQ6 + PC + COLQ + IL21R + DYRK1A, data = Top10T, family = "binomial")
summary(model)
oddsratio <- exp(coef(model)) 
write.csv(oddsratio,"OR_df9T.csv")
model_summary <- summary(model)
coeff_matrix <- model_summary$coefficients
coeff_df <- as.data.frame(coeff_matrix)
coeff_df$Term <- rownames(coeff_df)
rownames(coeff_df) <- NULL
coeff_df <- coeff_df[, c(ncol(coeff_df), 1:(ncol(coeff_df)-1))]
print(coeff_df)
write.csv(coeff_df,"coeff_df9T.csv")
predicted_probs <- predict(model, type = "response")
predicted_classes <- ifelse(predicted_probs > 0.5, 1, 0)
table(Predicted = predicted_classes, Actual = Top10T$Class)
predicted_probs <- predict(model, type = "response")
roc_curve <- roc(Top10T$Class, predicted_probs)
auc_value <- auc(roc_curve)
auc_value

# Fit logistic regression model for DNMT3B, TTLL11, TMEM192, WDR66, COQ6, PC, COLQ, IL21R, DYRK1A, PIH1D1
model <- glm(Class ~ DNMT3B + TTLL11 + TMEM192 + WDR66 + COQ6 + PC + COLQ + IL21R + DYRK1A + PIH1D1, data = Top10T, family = "binomial")
summary(model)
oddsratio <- exp(coef(model)) 
write.csv(oddsratio,"OR_df10T.csv")
model_summary <- summary(model)
coeff_matrix <- model_summary$coefficients
coeff_df <- as.data.frame(coeff_matrix)
coeff_df$Term <- rownames(coeff_df)
rownames(coeff_df) <- NULL
coeff_df <- coeff_df[, c(ncol(coeff_df), 1:(ncol(coeff_df)-1))]
print(coeff_df)
write.csv(coeff_df,"coeff_df10T.csv")
predicted_probs <- predict(model, type = "response")
predicted_classes <- ifelse(predicted_probs > 0.5, 1, 0)
table(Predicted = predicted_classes, Actual = Top10T$Class)
predicted_probs <- predict(model, type = "response")
roc_curve <- roc(Top10T$Class, predicted_probs)
auc_value <- auc(roc_curve)
auc_value

ggplot(Top10T, aes(x =Class, y =`DNMT3B`, fill=Class)) +
  geom_boxplot()+
  theme_classic(base_size=14)+
  theme(legend.position ="none")+
  labs(x = NULL, y = "DNMT3B Expression (log2 counts per million)")
ggsave("DNMT3B.png", height=5, width=2.5, dpi=320)

ggplot(Top10T, aes(x =Class, y =`DNMT3B`)) +
  geom_boxplot() +
  geom_point(aes(colour=Class)) + 
  theme_classic(base_size=18)+
  geom_text_repel(data = Top10T, aes(label = rownames(Top10T)))
ggsave("DNMT3B_IDs.png", dpi=320)

ggplot(Top10T, aes(x =Class, y =`TTLL11`, fill=Class)) +
  geom_boxplot()+
  theme_classic(base_size=14)+
  theme(legend.position ="none")+
  labs(x = NULL, y = "TTLL11 Expression (log2 counts per million)")
ggsave("TTLL11.png", height=5, width=2.5, dpi=320)

ggplot(Top10T, aes(x =Class, y =`TTLL11`)) +
  geom_boxplot() +
  geom_point(aes(colour=Class)) + 
  theme_classic(base_size=18)+
  geom_text_repel(data = Top10T, aes(label = rownames(Top10T)))
ggsave("TTLL11_IDs.png", dpi=320)

ggplot(Top10T, aes(x =Class, y =`TMEM192`, fill=Class)) +
  geom_boxplot()+
  theme_classic(base_size=14)+
  theme(legend.position ="none")+
  labs(x = NULL, y = "TMEM192 Expression (log2 counts per million)")
ggsave("TMEM192.png", height=5, width=2.5, dpi=320)

ggplot(Top10T, aes(x =Class, y =`TMEM192`)) +
  geom_boxplot() +
  geom_point(aes(colour=Class)) + 
  theme_classic(base_size=18)+
  geom_text_repel(data = Top10T, aes(label = rownames(Top10T)))
ggsave("TMEM192_IDs.png", dpi=320)

ggplot(Top10T, aes(x =Class, y =`WDR66`, fill=Class)) +
  geom_boxplot()+
  theme_classic(base_size=14)+
  theme(legend.position ="none")+
  labs(x = NULL, y = "WDR66 Expression (log2 counts per million)")
ggsave("WDR66.png", height=5, width=2.5, dpi=320)

ggplot(Top10T, aes(x =Class, y =`WDR66`)) +
  geom_boxplot() +
  geom_point(aes(colour=Class)) + 
  theme_classic(base_size=18)+
  geom_text_repel(data = Top10T, aes(label = rownames(Top10T)))
ggsave("WDR66_IDs.png", dpi=320)

ggplot(Top10T, aes(x =Class, y =`COQ6`, fill=Class)) +
  geom_boxplot()+
  theme_classic(base_size=14)+
  theme(legend.position ="none")+
  labs(x = NULL, y = "COQ6 Expression (log2 counts per million)")
ggsave("COQ6.png", height=5, width=2.5, dpi=320)

ggplot(Top10T, aes(x =Class, y =`COQ6`)) +
  geom_boxplot() +
  geom_point(aes(colour=Class)) + 
  theme_classic(base_size=18)+
  geom_text_repel(data = Top10T, aes(label = rownames(Top10T)))
ggsave("COQ6_IDs.png", dpi=320)

ggplot(Top10T, aes(x =Class, y =`IL21R`, fill=Class)) +
  geom_boxplot()+
  theme_classic(base_size=14)+
  theme(legend.position ="none")+
  labs(x = NULL, y = "IL21R Expression (log2 counts per million)")
ggsave("IL21R.png", height=5, width=2.5, dpi=320)

ggplot(Top10T, aes(x =Class, y =`IL21R`)) +
  geom_boxplot() +
  geom_point(aes(colour=Class)) + 
  theme_classic(base_size=18)+
  geom_text_repel(data = Top10T, aes(label = rownames(Top10T)))
ggsave("IL21R_IDs.png", dpi=320)

ggplot(Top10T, aes(x =Class, y =`DYRK1A`, fill=Class)) +
  geom_boxplot()+
  theme_classic(base_size=14)+
  theme(legend.position ="none")+
  labs(x = NULL, y = "DYRK1A Expression (log2 counts per million)")
ggsave("DYRK1A.png", height=5, width=2.5, dpi=320)

ggplot(Top10T, aes(x =Class, y =`DYRK1A`)) +
  geom_boxplot() +
  geom_point(aes(colour=Class)) + 
  theme_classic(base_size=18)+
  geom_text_repel(data = Top10T, aes(label = rownames(Top10T)))
ggsave("DYRK1A_IDs.png", dpi=320)

ggplot(Top10T, aes(x =Class, y =`PIH1D1`, fill=Class)) +
  geom_boxplot()+
  theme_classic(base_size=14)+
  theme(legend.position ="none")+
  labs(x = NULL, y = "PIH1D1 Expression (log2 counts per million)")
ggsave("PIH1D1.png", height=5, width=2.5, dpi=320)

ggplot(Top10T, aes(x =Class, y =`PIH1D1`)) +
  geom_boxplot() +
  geom_point(aes(colour=Class)) + 
  theme_classic(base_size=18)+
  geom_text_repel(data = Top10T, aes(label = rownames(Top10T)))
ggsave("PIH1D1_IDs.png", dpi=320)

# Fit logistic regression model for TTLL11
model <- glm(Class ~ TTLL11, data = Top10T, family = "binomial")
summary(model)
oddsratio <- exp(coef(model)) 
write.csv(oddsratio,"OR_dfTTLL11.csv")
model_summary <- summary(model)
coeff_matrix <- model_summary$coefficients
coeff_df <- as.data.frame(coeff_matrix)
coeff_df$Term <- rownames(coeff_df)
rownames(coeff_df) <- NULL
coeff_df <- coeff_df[, c(ncol(coeff_df), 1:(ncol(coeff_df)-1))]
print(coeff_df)
write.csv(coeff_df,"coeff_TTLL11.csv")
predicted_probs <- predict(model, type = "response")
predicted_classes <- ifelse(predicted_probs > 0.5, 1, 0)
table(Predicted = predicted_classes, Actual = Top10T$Class)
predicted_probs <- predict(model, type = "response")
roc_curve <- roc(Top10T$Class, predicted_probs)
auc_value <- auc(roc_curve)
auc_value

# Fit logistic regression model for TMEM192
model <- glm(Class ~ TMEM192, data = Top10T, family = "binomial")
summary(model)
oddsratio <- exp(coef(model)) 
write.csv(oddsratio,"OR_dfTMEM192.csv")
model_summary <- summary(model)
coeff_matrix <- model_summary$coefficients
coeff_df <- as.data.frame(coeff_matrix)
coeff_df$Term <- rownames(coeff_df)
rownames(coeff_df) <- NULL
coeff_df <- coeff_df[, c(ncol(coeff_df), 1:(ncol(coeff_df)-1))]
print(coeff_df)
write.csv(coeff_df,"coeff_TMEM192.csv")
predicted_probs <- predict(model, type = "response")
predicted_classes <- ifelse(predicted_probs > 0.5, 1, 0)
table(Predicted = predicted_classes, Actual = Top10T$Class)
predicted_probs <- predict(model, type = "response")
roc_curve <- roc(Top10T$Class, predicted_probs)
auc_value <- auc(roc_curve)
auc_value

# Fit logistic regression model for WDR66
model <- glm(Class ~ WDR66, data = Top10T, family = "binomial")
summary(model)
oddsratio <- exp(coef(model)) 
write.csv(oddsratio,"OR_dfWDR66.csv")
model_summary <- summary(model)
coeff_matrix <- model_summary$coefficients
coeff_df <- as.data.frame(coeff_matrix)
coeff_df$Term <- rownames(coeff_df)
rownames(coeff_df) <- NULL
coeff_df <- coeff_df[, c(ncol(coeff_df), 1:(ncol(coeff_df)-1))]
print(coeff_df)
write.csv(coeff_df,"coeff_WDR66.csv")
predicted_probs <- predict(model, type = "response")
predicted_classes <- ifelse(predicted_probs > 0.5, 1, 0)
table(Predicted = predicted_classes, Actual = Top10T$Class)
predicted_probs <- predict(model, type = "response")
roc_curve <- roc(Top10T$Class, predicted_probs)
auc_value <- auc(roc_curve)
auc_value

# Fit logistic regression model for COQ6
model <- glm(Class ~ COQ6, data = Top10T, family = "binomial")
summary(model)
oddsratio <- exp(coef(model)) 
write.csv(oddsratio,"OR_dfCOQ6.csv")
model_summary <- summary(model)
coeff_matrix <- model_summary$coefficients
coeff_df <- as.data.frame(coeff_matrix)
coeff_df$Term <- rownames(coeff_df)
rownames(coeff_df) <- NULL
coeff_df <- coeff_df[, c(ncol(coeff_df), 1:(ncol(coeff_df)-1))]
print(coeff_df)
write.csv(coeff_df,"coeff_COQ6.csv")
predicted_probs <- predict(model, type = "response")
predicted_classes <- ifelse(predicted_probs > 0.5, 1, 0)
table(Predicted = predicted_classes, Actual = Top10T$Class)
predicted_probs <- predict(model, type = "response")
roc_curve <- roc(Top10T$Class, predicted_probs)
auc_value <- auc(roc_curve)
auc_value

# Fit logistic regression model for PC
model <- glm(Class ~ PC, data = Top10T, family = "binomial")
summary(model)
oddsratio <- exp(coef(model)) 
write.csv(oddsratio,"OR_dfPC.csv")
model_summary <- summary(model)
coeff_matrix <- model_summary$coefficients
coeff_df <- as.data.frame(coeff_matrix)
coeff_df$Term <- rownames(coeff_df)
rownames(coeff_df) <- NULL
coeff_df <- coeff_df[, c(ncol(coeff_df), 1:(ncol(coeff_df)-1))]
print(coeff_df)
write.csv(coeff_df,"coeff_PC.csv")
predicted_probs <- predict(model, type = "response")
predicted_classes <- ifelse(predicted_probs > 0.5, 1, 0)
table(Predicted = predicted_classes, Actual = Top10T$Class)
predicted_probs <- predict(model, type = "response")
roc_curve <- roc(Top10T$Class, predicted_probs)
auc_value <- auc(roc_curve)
auc_value

# Fit logistic regression model for COLQ
model <- glm(Class ~ COLQ, data = Top10T, family = "binomial")
summary(model)
oddsratio <- exp(coef(model)) 
write.csv(oddsratio,"OR_dfCOLQ.csv") 
model_summary <- summary(model)
coeff_matrix <- model_summary$coefficients
coeff_df <- as.data.frame(coeff_matrix)
coeff_df$Term <- rownames(coeff_df)
rownames(coeff_df) <- NULL
coeff_df <- coeff_df[, c(ncol(coeff_df), 1:(ncol(coeff_df)-1))]
print(coeff_df)
write.csv(coeff_df,"coeff_COLQ.csv")
predicted_probs <- predict(model, type = "response")
predicted_classes <- ifelse(predicted_probs > 0.5, 1, 0)
table(Predicted = predicted_classes, Actual = Top10T$Class)
predicted_probs <- predict(model, type = "response")
roc_curve <- roc(Top10T$Class, predicted_probs)
auc_value <- auc(roc_curve)
auc_value

# Fit logistic regression model for IL21R
model <- glm(Class ~ IL21R, data = Top10T, family = "binomial")
summary(model)
oddsratio <- exp(coef(model)) 
write.csv(oddsratio,"OR_dfIL21R.csv") 
model_summary <- summary(model)
coeff_matrix <- model_summary$coefficients
coeff_df <- as.data.frame(coeff_matrix)
coeff_df$Term <- rownames(coeff_df)
rownames(coeff_df) <- NULL
coeff_df <- coeff_df[, c(ncol(coeff_df), 1:(ncol(coeff_df)-1))]
print(coeff_df)
write.csv(coeff_df,"coeff_IL21R.csv")
predicted_probs <- predict(model, type = "response")
predicted_classes <- ifelse(predicted_probs > 0.5, 1, 0)
table(Predicted = predicted_classes, Actual = Top10T$Class)
predicted_probs <- predict(model, type = "response")
roc_curve <- roc(Top10T$Class, predicted_probs)
auc_value <- auc(roc_curve)
auc_value

# Fit logistic regression model for DYRK1A
model <- glm(Class ~ DYRK1A, data = Top10T, family = "binomial")
summary(model)
oddsratio <- exp(coef(model)) 
write.csv(oddsratio,"OR_dfDYRK1A.csv") 
model_summary <- summary(model)
coeff_matrix <- model_summary$coefficients
coeff_df <- as.data.frame(coeff_matrix)
coeff_df$Term <- rownames(coeff_df)
rownames(coeff_df) <- NULL
coeff_df <- coeff_df[, c(ncol(coeff_df), 1:(ncol(coeff_df)-1))]
print(coeff_df)
write.csv(coeff_df,"coeff_DYRK1A.csv")
predicted_probs <- predict(model, type = "response")
predicted_classes <- ifelse(predicted_probs > 0.5, 1, 0)
table(Predicted = predicted_classes, Actual = Top10T$Class)
predicted_probs <- predict(model, type = "response")
roc_curve <- roc(Top10T$Class, predicted_probs)
auc_value <- auc(roc_curve)
auc_value

# Fit logistic regression model for PIH1D1
model <- glm(Class ~ PIH1D1, data = Top10T, family = "binomial")
summary(model)
oddsratio <- exp(coef(model)) 
write.csv(oddsratio,"OR_dfPIH1D1.csv") 
model_summary <- summary(model)
coeff_matrix <- model_summary$coefficients
coeff_df <- as.data.frame(coeff_matrix)
coeff_df$Term <- rownames(coeff_df)
rownames(coeff_df) <- NULL
coeff_df <- coeff_df[, c(ncol(coeff_df), 1:(ncol(coeff_df)-1))]
print(coeff_df)
write.csv(coeff_df,"coeff_PIH1D1.csv")
predicted_probs <- predict(model, type = "response")
predicted_classes <- ifelse(predicted_probs > 0.5, 1, 0)
table(Predicted = predicted_classes, Actual = Top10T$Class)
predicted_probs <- predict(model, type = "response")
roc_curve <- roc(Top10T$Class, predicted_probs)
auc_value <- auc(roc_curve)
auc_value

#Logistic regression - Methylome----
# Fit logistic regression model for 10 CpGs
model <- glm(Class ~ cg24387126 + cg27072813 + cg15626881 + cg22580353 + cg09877009 + cg08505473 + cg18984282 + cg15318627 + cg17421241 + cg25603742, data = Top10M, family = "binomial")
summary(model)
oddsratio <- exp(coef(model)) 
write.csv(oddsratio,"OR_10CpG.csv") 
model_summary <- summary(model)
coeff_matrix <- model_summary$coefficients
coeff_df <- as.data.frame(coeff_matrix)
coeff_df$Term <- rownames(coeff_df)
rownames(coeff_df) <- NULL
coeff_df <- coeff_df[, c(ncol(coeff_df), 1:(ncol(coeff_df)-1))]
print(coeff_df)
write.csv(coeff_df,"coeff_10CpG.csv")

predicted_probs <- predict(model, type = "response")
predicted_classes <- ifelse(predicted_probs > 0.5, 1, 0)

table(Predicted = predicted_classes, Actual = Top10M$Class)
predicted_probs <- predict(model, type = "response")
roc_curve <- roc(Top10M$Class, predicted_probs)
auc_value <- auc(roc_curve)
auc_value

# Fit logistic regression model for 1st CpG
model <- glm(Class ~ cg24387126, data = Top10M, family = "binomial")
summary(model)
oddsratio <- exp(coef(model)) 
write.csv(oddsratio,"OR_1stCpG.csv") 
model_summary <- summary(model)
coeff_matrix <- model_summary$coefficients
coeff_df <- as.data.frame(coeff_matrix)
coeff_df$Term <- rownames(coeff_df)
rownames(coeff_df) <- NULL
coeff_df <- coeff_df[, c(ncol(coeff_df), 1:(ncol(coeff_df)-1))]
print(coeff_df)
write.csv(coeff_df,"coeff_1stCpG.csv")

predicted_probs <- predict(model, type = "response")
predicted_classes <- ifelse(predicted_probs > 0.5, 1, 0)

table(Predicted = predicted_classes, Actual = Top10M$Class)
predicted_probs <- predict(model, type = "response")
roc_curve <- roc(Top10M$Class, predicted_probs)
auc_value <- auc(roc_curve)
auc_value

# Fit logistic regression model for 2nd CpG
model <- glm(Class ~ cg27072813, data = Top10M, family = "binomial")
summary(model)
oddsratio <- exp(coef(model)) 
write.csv(oddsratio,"OR_2ndCpG.csv") 
model_summary <- summary(model)
coeff_matrix <- model_summary$coefficients
coeff_df <- as.data.frame(coeff_matrix)
coeff_df$Term <- rownames(coeff_df)
rownames(coeff_df) <- NULL
coeff_df <- coeff_df[, c(ncol(coeff_df), 1:(ncol(coeff_df)-1))]
print(coeff_df)
write.csv(coeff_df,"coeff_2ndCpG.csv")

predicted_probs <- predict(model, type = "response")
predicted_classes <- ifelse(predicted_probs > 0.5, 1, 0)

table(Predicted = predicted_classes, Actual = Top10M$Class)
predicted_probs <- predict(model, type = "response")
roc_curve <- roc(Top10M$Class, predicted_probs)
auc_value <- auc(roc_curve)
auc_value

# Fit logistic regression model for 3rd CpG
model <- glm(Class ~ cg15626881, data = Top10M, family = "binomial")
summary(model)
oddsratio <- exp(coef(model)) 
write.csv(oddsratio,"OR_3rdCpG.csv") 
model_summary <- summary(model)
coeff_matrix <- model_summary$coefficients
coeff_df <- as.data.frame(coeff_matrix)
coeff_df$Term <- rownames(coeff_df)
rownames(coeff_df) <- NULL
coeff_df <- coeff_df[, c(ncol(coeff_df), 1:(ncol(coeff_df)-1))]
print(coeff_df)
write.csv(coeff_df,"coeff_3rdCpG.csv")

predicted_probs <- predict(model, type = "response")
predicted_classes <- ifelse(predicted_probs > 0.5, 1, 0)

table(Predicted = predicted_classes, Actual = Top10M$Class)
predicted_probs <- predict(model, type = "response")
roc_curve <- roc(Top10M$Class, predicted_probs)
auc_value <- auc(roc_curve)
auc_value

# Fit logistic regression model for 4th CpG
model <- glm(Class ~ cg22580353, data = Top10M, family = "binomial")
summary(model)
oddsratio <- exp(coef(model)) 
write.csv(oddsratio,"OR_4thCpG.csv") 
model_summary <- summary(model)
coeff_matrix <- model_summary$coefficients
coeff_df <- as.data.frame(coeff_matrix)
coeff_df$Term <- rownames(coeff_df)
rownames(coeff_df) <- NULL
coeff_df <- coeff_df[, c(ncol(coeff_df), 1:(ncol(coeff_df)-1))]
print(coeff_df)
write.csv(coeff_df,"coeff_4thCpG.csv")

predicted_probs <- predict(model, type = "response")
predicted_classes <- ifelse(predicted_probs > 0.5, 1, 0)

table(Predicted = predicted_classes, Actual = Top10M$Class)
predicted_probs <- predict(model, type = "response")
roc_curve <- roc(Top10M$Class, predicted_probs)
auc_value <- auc(roc_curve)
auc_value

# Fit logistic regression model for 5th CpG
model <- glm(Class ~ cg09877009, data = Top10M, family = "binomial")
summary(model)
oddsratio <- exp(coef(model)) 
write.csv(oddsratio,"OR_5thCpG.csv") 
model_summary <- summary(model)
coeff_matrix <- model_summary$coefficients
coeff_df <- as.data.frame(coeff_matrix)
coeff_df$Term <- rownames(coeff_df)
rownames(coeff_df) <- NULL
coeff_df <- coeff_df[, c(ncol(coeff_df), 1:(ncol(coeff_df)-1))]
print(coeff_df)
write.csv(coeff_df,"coeff_5thCpG.csv")

predicted_probs <- predict(model, type = "response")
predicted_classes <- ifelse(predicted_probs > 0.5, 1, 0)

table(Predicted = predicted_classes, Actual = Top10M$Class)
predicted_probs <- predict(model, type = "response")
roc_curve <- roc(Top10M$Class, predicted_probs)
auc_value <- auc(roc_curve)
auc_value

# Fit logistic regression model for 6th CpG
model <- glm(Class ~ cg08505473, data = Top10M, family = "binomial")
summary(model)
oddsratio <- exp(coef(model)) 
write.csv(oddsratio,"OR_6thCpG.csv") 
model_summary <- summary(model)
coeff_matrix <- model_summary$coefficients
coeff_df <- as.data.frame(coeff_matrix)
coeff_df$Term <- rownames(coeff_df)
rownames(coeff_df) <- NULL
coeff_df <- coeff_df[, c(ncol(coeff_df), 1:(ncol(coeff_df)-1))]
print(coeff_df)
write.csv(coeff_df,"coeff_6thCpG.csv")

predicted_probs <- predict(model, type = "response")
predicted_classes <- ifelse(predicted_probs > 0.5, 1, 0)

table(Predicted = predicted_classes, Actual = Top10M$Class)
predicted_probs <- predict(model, type = "response")
roc_curve <- roc(Top10M$Class, predicted_probs)
auc_value <- auc(roc_curve)
auc_value

# Fit logistic regression model for 7th CpG
model <- glm(Class ~ cg18984282, data = Top10M, family = "binomial")
summary(model)
oddsratio <- exp(coef(model)) 
write.csv(oddsratio,"OR_7thCpG.csv") 
model_summary <- summary(model)
coeff_matrix <- model_summary$coefficients
coeff_df <- as.data.frame(coeff_matrix)
coeff_df$Term <- rownames(coeff_df)
rownames(coeff_df) <- NULL
coeff_df <- coeff_df[, c(ncol(coeff_df), 1:(ncol(coeff_df)-1))]
print(coeff_df)
write.csv(coeff_df,"coeff_7thCpG.csv")
predicted_probs <- predict(model, type = "response")
predicted_classes <- ifelse(predicted_probs > 0.5, 1, 0)

table(Predicted = predicted_classes, Actual = Top10M$Class)
predicted_probs <- predict(model, type = "response")
roc_curve <- roc(Top10M$Class, predicted_probs)
auc_value <- auc(roc_curve)
auc_value

# Fit logistic regression model for 8th CpG
model <- glm(Class ~ cg15318627, data = Top10M, family = "binomial")
summary(model)
oddsratio <- exp(coef(model)) 
write.csv(oddsratio,"OR_8thCpG.csv") 
model_summary <- summary(model)
coeff_matrix <- model_summary$coefficients
coeff_df <- as.data.frame(coeff_matrix)
coeff_df$Term <- rownames(coeff_df)
rownames(coeff_df) <- NULL
coeff_df <- coeff_df[, c(ncol(coeff_df), 1:(ncol(coeff_df)-1))]
print(coeff_df)
write.csv(coeff_df,"coeff_8thCpG.csv")


predicted_probs <- predict(model, type = "response")
predicted_classes <- ifelse(predicted_probs > 0.5, 1, 0)

table(Predicted = predicted_classes, Actual = Top10M$Class)
predicted_probs <- predict(model, type = "response")
roc_curve <- roc(Top10M$Class, predicted_probs)
auc_value <- auc(roc_curve)
auc_value

# Fit logistic regression model for 9th CpG
model <- glm(Class ~ cg17421241, data = Top10M, family = "binomial")
summary(model)
oddsratio <- exp(coef(model)) 
write.csv(oddsratio,"OR_9thCpG.csv") 
model_summary <- summary(model)
coeff_matrix <- model_summary$coefficients
coeff_df <- as.data.frame(coeff_matrix)
coeff_df$Term <- rownames(coeff_df)
rownames(coeff_df) <- NULL
coeff_df <- coeff_df[, c(ncol(coeff_df), 1:(ncol(coeff_df)-1))]
print(coeff_df)
write.csv(coeff_df,"coeff_9thCpG.csv") 

predicted_probs <- predict(model, type = "response")
predicted_classes <- ifelse(predicted_probs > 0.5, 1, 0)

table(Predicted = predicted_classes, Actual = Top10M$Class)
predicted_probs <- predict(model, type = "response")
roc_curve <- roc(Top10M$Class, predicted_probs)
auc_value <- auc(roc_curve)
auc_value

# Fit logistic regression model for 10th CpG
model <- glm(Class ~ cg25603742, data = Top10M, family = "binomial")
summary(model)
oddsratio <- exp(coef(model)) 
write.csv(oddsratio,"OR_10thCpG.csv") 
model_summary <- summary(model)
coeff_matrix <- model_summary$coefficients
coeff_df <- as.data.frame(coeff_matrix)
coeff_df$Term <- rownames(coeff_df)
rownames(coeff_df) <- NULL
coeff_df <- coeff_df[, c(ncol(coeff_df), 1:(ncol(coeff_df)-1))]
print(coeff_df)
write.csv(coeff_df,"coeff_10thCpG.csv")


predicted_probs <- predict(model, type = "response")
predicted_classes <- ifelse(predicted_probs > 0.5, 1, 0)

table(Predicted = predicted_classes, Actual = Top10M$Class)
predicted_probs <- predict(model, type = "response")
roc_curve <- roc(Top10M$Class, predicted_probs)
auc_value <- auc(roc_curve)
auc_value

# Fit logistic regression model for CpG 1 & 2
model <- glm(Class ~ cg24387126 + cg27072813, data = Top10M, family = "binomial")
summary(model)
oddsratio <- exp(coef(model)) 
write.csv(oddsratio,"OR_2CpG.csv") 
model_summary <- summary(model)
coeff_matrix <- model_summary$coefficients
coeff_df <- as.data.frame(coeff_matrix)
coeff_df$Term <- rownames(coeff_df)
rownames(coeff_df) <- NULL
coeff_df <- coeff_df[, c(ncol(coeff_df), 1:(ncol(coeff_df)-1))]
print(coeff_df)
write.csv(coeff_df,"coeff_2CpG.csv")

predicted_probs <- predict(model, type = "response")
predicted_classes <- ifelse(predicted_probs > 0.5, 1, 0)

table(Predicted = predicted_classes, Actual = Top10M$Class)
predicted_probs <- predict(model, type = "response")
roc_curve <- roc(Top10M$Class, predicted_probs)
auc_value <- auc(roc_curve)
auc_value

# Fit logistic regression model for CpG 1, 2, 3
model <- glm(Class ~ cg24387126 + cg27072813 + cg15626881, data = Top10M, family = "binomial")
summary(model)
oddsratio <- exp(coef(model)) 
write.csv(oddsratio,"OR_3CpG.csv") 
model_summary <- summary(model)
coeff_matrix <- model_summary$coefficients
coeff_df <- as.data.frame(coeff_matrix)
coeff_df$Term <- rownames(coeff_df)
rownames(coeff_df) <- NULL
coeff_df <- coeff_df[, c(ncol(coeff_df), 1:(ncol(coeff_df)-1))]
print(coeff_df)
write.csv(coeff_df,"coeff_3CpG.csv")

predicted_probs <- predict(model, type = "response")
predicted_classes <- ifelse(predicted_probs > 0.5, 1, 0)

table(Predicted = predicted_classes, Actual = Top10M$Class)
predicted_probs <- predict(model, type = "response")
roc_curve <- roc(Top10M$Class, predicted_probs)
auc_value <- auc(roc_curve)
auc_value

# Fit logistic regression model for CpG 1-4
model <- glm(Class ~ cg24387126 + cg27072813 + cg15626881 + cg22580353, data = Top10M, family = "binomial")
summary(model)
oddsratio <- exp(coef(model)) 
write.csv(oddsratio,"OR_4CpG.csv") 
model_summary <- summary(model)
coeff_matrix <- model_summary$coefficients
coeff_df <- as.data.frame(coeff_matrix)
coeff_df$Term <- rownames(coeff_df)
rownames(coeff_df) <- NULL
coeff_df <- coeff_df[, c(ncol(coeff_df), 1:(ncol(coeff_df)-1))]
print(coeff_df)
write.csv(coeff_df,"coeff_4CpG.csv")

predicted_probs <- predict(model, type = "response")
predicted_classes <- ifelse(predicted_probs > 0.5, 1, 0)

table(Predicted = predicted_classes, Actual = Top10M$Class)
predicted_probs <- predict(model, type = "response")
roc_curve <- roc(Top10M$Class, predicted_probs)
auc_value <- auc(roc_curve)
auc_value

# Fit logistic regression model for CpG 1-5
model <- glm(Class ~ cg24387126 + cg27072813 + cg15626881 + cg22580353 + cg09877009, data = Top10M, family = "binomial")
summary(model)
oddsratio <- exp(coef(model)) 
write.csv(oddsratio,"OR_5CpG.csv") 
model_summary <- summary(model)
coeff_matrix <- model_summary$coefficients
coeff_df <- as.data.frame(coeff_matrix)
coeff_df$Term <- rownames(coeff_df)
rownames(coeff_df) <- NULL
coeff_df <- coeff_df[, c(ncol(coeff_df), 1:(ncol(coeff_df)-1))]
print(coeff_df)
write.csv(coeff_df,"coeff_5CpG.csv")

predicted_probs <- predict(model, type = "response")
predicted_classes <- ifelse(predicted_probs > 0.5, 1, 0)

table(Predicted = predicted_classes, Actual = Top10M$Class)
predicted_probs <- predict(model, type = "response")
roc_curve <- roc(Top10M$Class, predicted_probs)
auc_value <- auc(roc_curve)
auc_value

# Fit logistic regression model for CpG 1-6
model <- glm(Class ~ cg24387126 + cg27072813 + cg15626881 + cg22580353 + cg09877009 + cg08505473, data = Top10M, family = "binomial")
summary(model)
oddsratio <- exp(coef(model)) 
write.csv(oddsratio,"OR_6CpG.csv") 
model_summary <- summary(model)
coeff_matrix <- model_summary$coefficients
coeff_df <- as.data.frame(coeff_matrix)
coeff_df$Term <- rownames(coeff_df)
rownames(coeff_df) <- NULL
coeff_df <- coeff_df[, c(ncol(coeff_df), 1:(ncol(coeff_df)-1))]
print(coeff_df)
write.csv(coeff_df,"coeff_6CpG.csv") 

predicted_probs <- predict(model, type = "response")
predicted_classes <- ifelse(predicted_probs > 0.5, 1, 0)

table(Predicted = predicted_classes, Actual = Top10M$Class)
predicted_probs <- predict(model, type = "response")
roc_curve <- roc(Top10M$Class, predicted_probs)
auc_value <- auc(roc_curve)
auc_value

# Fit logistic regression model for CpG 1-7
model <- glm(Class ~ cg24387126 + cg27072813 + cg15626881 + cg22580353 + cg09877009 + cg08505473 + cg18984282, data = Top10M, family = "binomial")
summary(model)
oddsratio <- exp(coef(model)) 
write.csv(oddsratio,"OR_7CpG.csv") 
model_summary <- summary(model)
coeff_matrix <- model_summary$coefficients
coeff_df <- as.data.frame(coeff_matrix)
coeff_df$Term <- rownames(coeff_df)
rownames(coeff_df) <- NULL
coeff_df <- coeff_df[, c(ncol(coeff_df), 1:(ncol(coeff_df)-1))]
print(coeff_df)
write.csv(coeff_df,"coeff_7CpG.csv")

predicted_probs <- predict(model, type = "response")
predicted_classes <- ifelse(predicted_probs > 0.5, 1, 0)

table(Predicted = predicted_classes, Actual = Top10M$Class)
predicted_probs <- predict(model, type = "response")
roc_curve <- roc(Top10M$Class, predicted_probs)
auc_value <- auc(roc_curve)
auc_value

# Fit logistic regression model for CpG 1-8
model <- glm(Class ~ cg24387126 + cg27072813 + cg15626881 + cg22580353 + cg09877009 + cg08505473 + cg18984282 + cg15318627, data = Top10M, family = "binomial")
summary(model)
oddsratio <- exp(coef(model)) 
write.csv(oddsratio,"OR_8CpG.csv") 
model_summary <- summary(model)
coeff_matrix <- model_summary$coefficients
coeff_df <- as.data.frame(coeff_matrix)
coeff_df$Term <- rownames(coeff_df)
rownames(coeff_df) <- NULL
coeff_df <- coeff_df[, c(ncol(coeff_df), 1:(ncol(coeff_df)-1))]
print(coeff_df)
write.csv(coeff_df,"coeff_8CpG.csv")

predicted_probs <- predict(model, type = "response")
predicted_classes <- ifelse(predicted_probs > 0.5, 1, 0)

table(Predicted = predicted_classes, Actual = Top10M$Class)
predicted_probs <- predict(model, type = "response")
roc_curve <- roc(Top10M$Class, predicted_probs)
auc_value <- auc(roc_curve)
auc_value

# Fit logistic regression model for CpG 1-9
model <- glm(Class ~ cg24387126 + cg27072813 + cg15626881 + cg22580353 + cg09877009 + cg08505473 + cg18984282 + cg15318627 + cg17421241, data = Top10M, family = "binomial")
summary(model)
oddsratio <- exp(coef(model)) 
write.csv(oddsratio,"OR_9CpG.csv") 
model_summary <- summary(model)
coeff_matrix <- model_summary$coefficients
coeff_df <- as.data.frame(coeff_matrix)
coeff_df$Term <- rownames(coeff_df)
rownames(coeff_df) <- NULL
coeff_df <- coeff_df[, c(ncol(coeff_df), 1:(ncol(coeff_df)-1))]
print(coeff_df)
write.csv(coeff_df,"coeff_9CpG.csv")

predicted_probs <- predict(model, type = "response")
predicted_classes <- ifelse(predicted_probs > 0.5, 1, 0)

table(Predicted = predicted_classes, Actual = Top10M$Class)
predicted_probs <- predict(model, type = "response")
roc_curve <- roc(Top10M$Class, predicted_probs)
auc_value <- auc(roc_curve)
auc_value

#Boxplots----

ggplot(Top10M, aes(x =Class, y =`cg24387126`, fill=Class)) +
  geom_boxplot()+
  theme_classic(base_size=14)+
  theme(legend.position ="none")+
  labs(x = NULL, y = "cg24387126 (Normalised Beta Values)")
ggsave("cg24387126.png", height=5, width=2.5, dpi=320)

ggplot(Top10M, aes(x =Class, y =`cg24387126`)) +
  geom_boxplot() +
  geom_point(aes(colour=Class)) + 
  theme_classic(base_size=18)+
  geom_text_repel(data = Top10M, aes(label = rownames(Top10M)))
ggsave("cg24387126_IDs.png", dpi=320)

ggplot(Top10M, aes(x =Class, y =`cg27072813`, fill=Class)) +
  geom_boxplot()+
  theme_classic(base_size=14)+
  theme(legend.position ="none")+
  labs(x = NULL, y = "cg24387126 (Normalised Beta Values)")
ggsave("cg27072813.png", height=5, width=2.5, dpi=320)

ggplot(Top10M, aes(x =Class, y =`cg27072813`)) +
  geom_boxplot() +
  geom_point(aes(colour=Class)) + 
  theme_classic(base_size=18)+
  geom_text_repel(data = Top10M, aes(label = rownames(Top10M)))
ggsave("cg27072813_IDs.png", dpi=320)

ggplot(Top10M, aes(x =Class, y =`cg15626881`, fill=Class)) +
  geom_boxplot()+
  theme_classic(base_size=14)+
  theme(legend.position ="none")+
  labs(x = NULL, y = "cg15626881 (Normalised Beta Values)")
ggsave("cg15626881.png", height=5, width=2.5, dpi=320)

ggplot(Top10M, aes(x =Class, y =`cg15626881`)) +
  geom_boxplot() +
  geom_point(aes(colour=Class)) + 
  theme_classic(base_size=18)+
  geom_text_repel(data = Top10M, aes(label = rownames(Top10M)))
ggsave("cg15626881_IDs.png", dpi=320)

ggplot(Top10M, aes(x =Class, y =`cg22580353`, fill=Class)) +
  geom_boxplot()+
  theme_classic(base_size=14)+
  theme(legend.position ="none")+
  labs(x = NULL, y = "cg22580353 (Normalised Beta Values)")
ggsave("cg22580353.png", height=5, width=2.5, dpi=320)

ggplot(Top10M, aes(x =Class, y =`cg22580353`)) +
  geom_boxplot() +
  geom_point(aes(colour=Class)) + 
  theme_classic(base_size=18)+
  geom_text_repel(data = Top10M, aes(label = rownames(Top10M)))
ggsave("cg22580353_IDs.png", dpi=320)

ggplot(Top10M, aes(x =Class, y =`cg09877009`, fill=Class)) +
  geom_boxplot()+
  theme_classic(base_size=14)+
  theme(legend.position ="none")+
  labs(x = NULL, y = "cg09877009 (Normalised Beta Values)")
ggsave("cg09877009.png", height=5, width=2.5, dpi=320)

ggplot(Top10M, aes(x =Class, y =`cg09877009`)) +
  geom_boxplot() +
  geom_point(aes(colour=Class)) + 
  theme_classic(base_size=18)+
  geom_text_repel(data = Top10M, aes(label = rownames(Top10M)))
ggsave("cg09877009_IDs.png", dpi=320)

ggplot(Top10M, aes(x =Class, y =`cg08505473`, fill=Class)) +
  geom_boxplot()+
  theme_classic(base_size=14)+
  theme(legend.position ="none")+
  labs(x = NULL, y = "cg08505473 (Normalised Beta Values)")
ggsave("cg08505473.png", height=5, width=2.5, dpi=320)

ggplot(Top10M, aes(x =Class, y =`cg08505473`)) +
  geom_boxplot() +
  geom_point(aes(colour=Class)) + 
  theme_classic(base_size=18)+
  geom_text_repel(data = Top10M, aes(label = rownames(Top10M)))
ggsave("cg08505473_IDs.png", dpi=320)

ggplot(Top10M, aes(x =Class, y =`cg18984282`, fill=Class)) +
  geom_boxplot()+
  theme_classic(base_size=14)+
  theme(legend.position ="none")+
  labs(x = NULL, y = "cg18984282 (Normalised Beta Values)")
ggsave("cg18984282.png", height=5, width=2.5, dpi=320)

ggplot(Top10M, aes(x =Class, y =`cg18984282`)) +
  geom_boxplot() +
  geom_point(aes(colour=Class)) + 
  theme_classic(base_size=18)+
  geom_text_repel(data = Top10M, aes(label = rownames(Top10M)))
ggsave("cg18984282_IDs.png", dpi=320)

ggplot(Top10M, aes(x =Class, y =`cg15318627`, fill=Class)) +
  geom_boxplot()+
  theme_classic(base_size=14)+
  theme(legend.position ="none")+
  labs(x = NULL, y = "cg15318627 (Normalised Beta Values)")
ggsave("cg15318627.png", height=5, width=2.5, dpi=320)

ggplot(Top10M, aes(x =Class, y =`cg15318627`)) +
  geom_boxplot() +
  geom_point(aes(colour=Class)) + 
  theme_classic(base_size=18)+
  geom_text_repel(data = Top10M, aes(label = rownames(Top10M)))
ggsave("cg15318627_IDs.png", dpi=320)

ggplot(Top10M, aes(x =Class, y =`cg17421241`, fill=Class)) +
  geom_boxplot()+
  theme_classic(base_size=14)+
  theme(legend.position ="none")+
  labs(x = NULL, y = "cg17421241 (Normalised Beta Values)")
ggsave("cg17421241.png", height=5, width=2.5, dpi=320)

ggplot(Top10M, aes(x =Class, y =`cg17421241`)) +
  geom_boxplot() +
  geom_point(aes(colour=Class)) + 
  theme_classic(base_size=18)+
  geom_text_repel(data = Top10M, aes(label = rownames(Top10M)))
ggsave("cg17421241_IDs.png", dpi=320)

ggplot(Top10M, aes(x =Class, y =`cg25603742`, fill=Class)) +
  geom_boxplot()+
  theme_classic(base_size=14)+
  theme(legend.position ="none")+
  labs(x = NULL, y = "cg25603742 (Normalised Beta Values)")
ggsave("cg25603742.png", height=5, width=2.5, dpi=320)

ggplot(Top10M, aes(x =Class, y =`cg25603742`)) +
  geom_boxplot() +
  geom_point(aes(colour=Class)) + 
  theme_classic(base_size=18)+
  geom_text_repel(data = Top10M, aes(label = rownames(Top10M)))
ggsave("cg25603742_IDs.png", dpi=320)

