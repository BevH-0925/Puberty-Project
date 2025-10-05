setwd("~/Documents/Code/RMCH_cross_validation")

library(caret)
library(randomForest)
library(kernelshap)
library(shapviz)
library(pROC)
library(tidyverse)
library(ggrepel)
library(Boruta)

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

data<-cbind(Top500,Class)

set.seed(123)
# Run Boruta 
boruta_output <- Boruta(Class ~., data = data, doTrace = 2)

# Check Confirmed Attributes
final_features <- getSelectedAttributes(boruta_output, withTentative = FALSE)

# Subset data to include only the confirmed features + target
data_selected <- data[, c(final_features, "Class")]

# #Define stratified k-fold cross-validation
ctrl <- trainControl(method = "cv", number = 5,
                     classProbs = TRUE,
                     summaryFunction = twoClassSummary,
                     savePredictions = "final")

#Train Random Forest model with stratified CV
set.seed(123)
rf_model <- train(Class ~., data = data_selected,
                  method = "rf",
                  trControl = ctrl,
                  metric = "ROC",
                  tuneLength = 5)

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
oob_auc
final_features
write.csv (final_features, "Borutagenes_RMCH_glm_MF.csv")

resampled <- rf_model[["resampledCM"]]
write.csv(resampled, "FoldperformanceRMCHT_MFGLM.csv")
pred <- rf_model[["pred"]]
write.csv(pred, "PredRMCHT_MFGLM.csv")

#Shapley Values----
set.seed(123)
s <- kernelshap(rf_model, data_selected[-10], bg_X=data_selected, type="prob")
sv <- shapviz(s)
png(filename = "SHAPPlot_T_BorutaImp_glm_MF.png", width = 6 * 320, height = 4 * 320, res = 320)
sv_importance(sv, kind = "bee")
dev.off()

