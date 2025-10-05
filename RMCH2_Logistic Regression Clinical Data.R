setwd("~/Documents/Code/RMCH Clinical Data")
library(pROC)
library(tidyverse)

set.seed(123)

clindata <- read.csv("~/Documents/Code/RMCH Clinical Data/Data/Clindata_all.csv")
clindata$Group <- as.factor(clindata$Group)

# Fit logistic regression model for LH 0min
clindata$Group <- relevel(clindata$Group, ref = "Pubertal")
model <- glm(Group ~ LH0min, data = clindata, family = "binomial")
summary(model)
oddsratio <- exp(coef(model))
oddsratio
write.csv(oddsratio,"OR_LH0min.csv")

model_summary <- summary(model)
coeff_matrix <- model_summary$coefficients
coeff_df <- as.data.frame(coeff_matrix)
coeff_df$Term <- rownames(coeff_df)
rownames(coeff_df) <- NULL
coeff_df <- coeff_df[, c(ncol(coeff_df), 1:(ncol(coeff_df)-1))]
write.csv(coeff_df,"coeff_LH0min.csv")

predicted_probs <- predict(model, type = "response")
predicted_classes <- ifelse(predicted_probs > 0.5, 1, 0)

table(Predicted = predicted_classes, Actual = clindata$Group)
predicted_probs <- predict(model, type = "response")
roc_curveLH0 <- roc(clindata$Group, predicted_probs)
auc_valueLH0 <- auc(roc_curveLH0)
auc_valueLH0

best_coordsLH0 <- coords(roc_curveLH0, "best", ret = c("threshold", "sensitivity", "specificity"), transpose = FALSE)
best_coordsLH0

#To translate back to cut-off in original variable
logit <- function(p) log(p / (1 - p))
cutoff_value <- (logit(best_coordsLH0[1]) - coef(model)[1]) / coef(model)[2]
cutoff_value

ggroc(roc_curveLH0, colour = 'steelblue', linewidth = 1, legacy.axes = TRUE) +
  ggtitle(paste0('ROC Curve ', '(AUC = ', round(auc_valueLH0,2), ')')) +
  geom_abline(slope=1,intercept=0,lwd=1,lty=2,col="gray")+
  theme_classic(base_size = 16)+
  labs(x="False Positive Rate", y="True Positive Rate")+
  annotate("point", x = 1 - best_coordsLH0$specificity, y = best_coordsLH0$sensitivity, colour = "red", size = 3) +
  annotate("text", x = 1 - best_coordsLH0$specificity, y = best_coordsLH0$sensitivity, 
           label = paste("Cut-off =", round(cutoff_value, 1), "IU/L"), 
           vjust = -2, hjust = 0.2, colour = "red")
ggsave("ROC_LH0min.png", dpi=320)

# Fit logistic regression model for FSH 0min
clindata$Group <- relevel(clindata$Group, ref = "Pubertal")
model <- glm(Group ~ FSH0min, data = clindata, family = "binomial")
summary(model)
oddsratio <- exp(coef(model))
oddsratio
write.csv(oddsratio,"OR_FSH0min.csv")

model_summary <- summary(model)
coeff_matrix <- model_summary$coefficients
coeff_df <- as.data.frame(coeff_matrix)
coeff_df$Term <- rownames(coeff_df)
rownames(coeff_df) <- NULL
coeff_df <- coeff_df[, c(ncol(coeff_df), 1:(ncol(coeff_df)-1))]
write.csv(coeff_df,"coeff_FSH0min.csv")

predicted_probs <- predict(model, type = "response")
predicted_classes <- ifelse(predicted_probs > 0.5, 1, 0)

table(Predicted = predicted_classes, Actual = clindata$Group)
predicted_probs <- predict(model, type = "response")
roc_curveFSH0 <- roc(clindata$Group, predicted_probs)
auc_valueFSH0 <- auc(roc_curveFSH0)
auc_valueFSH0

best_coordsFSH0 <- coords(roc_curveFSH0, "best", ret = c("threshold", "sensitivity", "specificity"), transpose = FALSE)
best_coordsFSH0

#To translate back to cut-off in original variable
logit <- function(p) log(p / (1 - p))
cutoff_value <- (logit(best_coordsFSH0[1]) - coef(model)[1]) / coef(model)[2]
cutoff_value

ggroc(roc_curveFSH0, colour = 'steelblue', linewidth = 1, legacy.axes = TRUE) +
  ggtitle(paste0('ROC Curve ', '(AUC = ', round(auc_valueFSH0,2), ')')) +
  geom_abline(slope=1,intercept=0,lwd=1,lty=2,col="gray")+
  theme_classic(base_size = 16)+
  labs(x="False Positive Rate", y="True Positive Rate")+
  annotate("point", x = 1 - best_coordsFSH0$specificity, y = best_coordsFSH0$sensitivity, colour = "red", size = 3) +
  annotate("text", x = 1 - best_coordsFSH0$specificity, y = best_coordsFSH0$sensitivity, 
           label = paste("Cut-off =", round(cutoff_value, 1), "IU/L"), 
           vjust = -1, hjust = 0.6, colour = "red")
ggsave("ROC_FSH0min.png", dpi=320)

# Fit logistic regression model for Oestradiol
#Subset of clinical data with oestradiol concs
clindataoest <- clindata %>% drop_na(OestradiolConc)
clindataoest$Group <- relevel(clindataoest$Group, ref = "Pubertal")

model <- glm(Group ~ OestradiolConc, data = clindataoest, family = "binomial")
summary(model)
oddsratio <- exp(coef(model))
oddsratio
write.csv(oddsratio,"OR_Oest.csv")

model_summary <- summary(model)
coeff_matrix <- model_summary$coefficients
coeff_df <- as.data.frame(coeff_matrix)
coeff_df$Term <- rownames(coeff_df)
rownames(coeff_df) <- NULL
coeff_df <- coeff_df[, c(ncol(coeff_df), 1:(ncol(coeff_df)-1))]
write.csv(coeff_df,"coeff_Oest.csv")

predicted_probs <- predict(model, type = "response")
predicted_classes <- ifelse(predicted_probs > 0.5, 1, 0)

table(Predicted = predicted_classes, Actual = clindataoest$Group)
predicted_probs <- predict(model, type = "response")
roc_curveOest <- roc(clindataoest$Group, predicted_probs)
auc_valueOest<- auc(roc_curveOest)
auc_valueOest

best_coordsOest <- coords(roc_curveOest, "best", ret = c("threshold", "sensitivity", "specificity"), transpose = FALSE)
best_coordsOest

#To translate back to cut-off in original variable
logit <- function(p) log(p / (1 - p))
cutoff_value <- (logit(best_coordsOest[1]) - coef(model)[1]) / coef(model)[2]
cutoff_value

ggroc(roc_curveOest, colour = 'steelblue', linewidth = 1, legacy.axes = TRUE) +
  ggtitle(paste0('ROC Curve ', '(AUC = ', round(auc_valueOest,2), ')')) +
  geom_abline(slope=1,intercept=0,lwd=1,lty=2,col="gray")+
  theme_classic(base_size = 16)+
  labs(x="False Positive Rate", y="True Positive Rate")+
  annotate("point", x = 1 - best_coordsOest$specificity, y = best_coordsOest$sensitivity, colour = "red", size = 3) +
  annotate("text", x = 1 - best_coordsOest$specificity, y = best_coordsOest$sensitivity, 
           label = paste("Cut-off =", round(cutoff_value, 0), "pmol/L"), 
           vjust = -2, hjust = 0.1, colour = "red")
ggsave("ROC_Oest.png", dpi=320)

# Fit logistic regression model for testosterone
#Subset of clinical data with testosteroneconcs
clindatatesto <- clindata %>% drop_na(TestosteroneConc)
clindatatesto$Group <- relevel(clindatatesto$Group, ref = "Pubertal")
model <- glm(Group ~TestosteroneConc,  data = clindatatesto, family = "binomial")
summary(model)
oddsratio <- exp(coef(model))
oddsratio
write.csv(oddsratio,"OR_Testo.csv")

model_summary <- summary(model)
coeff_matrix <- model_summary$coefficients
coeff_df <- as.data.frame(coeff_matrix)
coeff_df$Term <- rownames(coeff_df)
rownames(coeff_df) <- NULL
coeff_df <- coeff_df[, c(ncol(coeff_df), 1:(ncol(coeff_df)-1))]
write.csv(coeff_df,"coeff_Oest.csv")

predicted_probs <- predict(model, type = "response")
predicted_classes <- ifelse(predicted_probs > 0.5, 1, 0)

table(Predicted = predicted_classes, Actual = clindatatesto$Group)
predicted_probs <- predict(model, type = "response")
roc_curveTesto<- roc(clindatatesto$Group, predicted_probs)
auc_valueTesto<- auc(roc_curveTesto)
auc_valueTesto

best_coordsTesto<- coords(roc_curveTesto, "best", ret = c("threshold", "sensitivity", "specificity"), transpose = FALSE)
best_coordsTesto

#To translate back to cut-off in original variable
logit <- function(p) log(p / (1 - p))
cutoff_value <- (logit(best_coordsTesto[1]) - coef(model)[1]) / coef(model)[2]
cutoff_value

ggroc(roc_curveTesto, colour = 'steelblue', linewidth = 1, legacy.axes = TRUE) +
  ggtitle(paste0('ROC Curve ', '(AUC = ', round(auc_valueTesto,2), ')')) +
  geom_abline(slope=1,intercept=0,lwd=1,lty=2,col="gray")+
  theme_classic(base_size = 16)+
  labs(x="False Positive Rate", y="True Positive Rate")+
  annotate("point", x = 1 - best_coordsTesto$specificity, y = best_coordsTesto$sensitivity, colour = "red", size = 3) +
  annotate("text", x = 1 - best_coordsTesto$specificity, y = best_coordsTesto$sensitivity, 
           label = paste("Cut-off =", round(cutoff_value, 1), "nmol/L"), 
           vjust = -2, hjust = 0.1, colour = "red")
ggsave("ROC_Testo.png", dpi=320)

# Fit logistic regression model for basal LH/FSH ratio
#Subset of clinical data with testosteroneconcs
clindata$Group <- relevel(clindata$Group, ref = "Pubertal")
model <- glm(Group ~BasalLHFSHratio,  data = clindata, family = "binomial")
summary(model)
oddsratio <- exp(coef(model))
oddsratio
write.csv(oddsratio,"OR_basalLHFSHratio.csv")

model_summary <- summary(model)
coeff_matrix <- model_summary$coefficients
coeff_df <- as.data.frame(coeff_matrix)
coeff_df$Term <- rownames(coeff_df)
rownames(coeff_df) <- NULL
coeff_df <- coeff_df[, c(ncol(coeff_df), 1:(ncol(coeff_df)-1))]
write.csv(coeff_df,"coeff_basalLHFSHratio.csv")

predicted_probs <- predict(model, type = "response")
predicted_classes <- ifelse(predicted_probs > 0.5, 1, 0)

table(Predicted = predicted_classes, Actual = clindata$Group)
predicted_probs <- predict(model, type = "response")
roc_curvebasalLHFSHratio<- roc(clindata$Group, predicted_probs)
auc_valuebasalLHFSHratio<- auc(roc_curvebasalLHFSHratio)
auc_valuebasalLHFSHratio

best_coordsbasalLHFSHratio<- coords(roc_curvebasalLHFSHratio, "best", ret = c("threshold", "sensitivity", "specificity"), transpose = FALSE)
best_coordsbasalLHFSHratio

#To translate back to cut-off in original variable
logit <- function(p) log(p / (1 - p))
cutoff_value <- (logit(best_coordsbasalLHFSHratio[1]) - coef(model)[1]) / coef(model)[2]
cutoff_value

ggroc(roc_curvebasalLHFSHratio, colour = 'steelblue', linewidth = 1, legacy.axes = TRUE) +
  ggtitle(paste0('ROC Curve ', '(AUC = ', round(auc_valuebasalLHFSHratio,2), ')')) +
  geom_abline(slope=1,intercept=0,lwd=1,lty=2,col="gray")+
  theme_classic(base_size = 16)+
  labs(x="False Positive Rate", y="True Positive Rate")+
  annotate("point", x = 1 - best_coordsbasalLHFSHratio$specificity, y = best_coordsbasalLHFSHratio$sensitivity, colour = "red", size = 3) +
  annotate("text", x = 1 - best_coordsbasalLHFSHratio$specificity, y = best_coordsbasalLHFSHratio$sensitivity, 
           label = paste("Cut-off =", round(cutoff_value, 1)), 
           vjust = -2, hjust = 0.1, colour = "red")
ggsave("ROC_basalLHFSHratio.png", dpi=320)

# Combined
roc_list <- list("Testosterone (AUC=0.97)"=roc_curveTesto, "Oestradiol (AUC=0.97)" = roc_curveOest, "Basal LH (AUC=0.95)" = roc_curveLH0, "Basal FSH (AUC=0.86)" = roc_curveFSH0, "Basal LH:FSH ratio (AUC=0.79)"=roc_curvebasalLHFSHratio)
ggroc(roc_list, aes = c("colour"), legacy.axes = TRUE)+
  geom_abline(slope=1,intercept=0,lwd=1,lty=2,col="gray")+
  theme_classic(base_size = 16)+
  labs(x="False Positive Rate", y="True Positive Rate", colour="Hormone")
ggsave("ROC_combined.png", dpi=320)

# Fit logistic regression model for basal LH & oestradiol
clindataoest$Group <- relevel(clindataoest$Group, ref = "Pubertal")
model <- glm(Group ~LH0min+OestradiolConc,  data = clindataoest, family = "binomial")
summary(model)
oddsratio <- exp(coef(model))
oddsratio
write.csv(oddsratio,"OR_LHoest.csv")

model_summary <- summary(model)
coeff_matrix <- model_summary$coefficients
coeff_df <- as.data.frame(coeff_matrix)
coeff_df$Term <- rownames(coeff_df)
rownames(coeff_df) <- NULL
coeff_df <- coeff_df[, c(ncol(coeff_df), 1:(ncol(coeff_df)-1))]
write.csv(coeff_df,"coeff_LHoest.csv")

predicted_probs <- predict(model, type = "response")
predicted_classes <- ifelse(predicted_probs > 0.5, 1, 0)

table(Predicted = predicted_classes, Actual = clindataoest$Group)
predicted_probs <- predict(model, type = "response")
roc_curveLHoest<- roc(clindataoest$Group, predicted_probs)
auc_valueLHoest<- auc(roc_curveLHoest)
auc_valueLHoest

best_coordsLHoest<- coords(roc_curveLHoest, "best", ret = c("threshold", "sensitivity", "specificity"), transpose = FALSE)
best_coordsLHoest

ggroc(roc_curveLHoest, colour = 'steelblue', linewidth = 1, legacy.axes = TRUE) +
  ggtitle(paste0('ROC Curve ', '(AUC = ', round(auc_valueLHoest,2), ')')) +
  geom_abline(slope=1,intercept=0,lwd=1,lty=2,col="gray")+
  theme_classic(base_size = 16)+
  labs(x="False Positive Rate", y="True Positive Rate")+
  annotate("point", x = 1 - best_coordsLHoest$specificity, y = best_coordsLHoest$sensitivity, colour = "red", size = 3)
ggsave("ROC_LHoest.png", dpi=320)

# Fit logistic regression model for basal LH & testo
clindatatesto$Group <- relevel(clindatatesto$Group, ref = "Pubertal")
model <- glm(Group ~LH0min+TestosteroneConc,  data = clindatatesto, family = "binomial")
summary(model)
oddsratio <- exp(coef(model))
oddsratio
write.csv(oddsratio,"OR_LHtesto.csv")

model_summary <- summary(model)
coeff_matrix <- model_summary$coefficients
coeff_df <- as.data.frame(coeff_matrix)
coeff_df$Term <- rownames(coeff_df)
rownames(coeff_df) <- NULL
coeff_df <- coeff_df[, c(ncol(coeff_df), 1:(ncol(coeff_df)-1))]
write.csv(coeff_df,"coeff_LHtesto.csv")

predicted_probs <- predict(model, type = "response")
predicted_classes <- ifelse(predicted_probs > 0.5, 1, 0)

table(Predicted = predicted_classes, Actual = clindatatesto$Group)
predicted_probs <- predict(model, type = "response")
roc_curveLHtesto<- roc(clindatatesto$Group, predicted_probs)
auc_valueLHtesto<- auc(roc_curveLHtesto)
auc_valueLHtesto

best_coordsLHtesto<- coords(roc_curveLHtesto, "best", ret = c("threshold", "sensitivity", "specificity"), transpose = FALSE)
best_coordsLHtesto

ggroc(roc_curveLHtesto, colour = 'steelblue', linewidth = 1, legacy.axes = TRUE) +
  ggtitle(paste0('ROC Curve ', '(AUC = ', round(auc_valueLHtesto,2), ')')) +
  geom_abline(slope=1,intercept=0,lwd=1,lty=2,col="gray")+
  theme_classic(base_size = 16)+
  labs(x="False Positive Rate", y="True Positive Rate")+
  annotate("point", x = 1 - best_coordsLHtesto$specificity, y = best_coordsLHtesto$sensitivity, colour = "red", size = 3)
ggsave("ROC_LHtesto.png", dpi=320)

