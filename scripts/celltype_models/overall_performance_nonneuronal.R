# Setup ####
## libraries
library(randomForest)
library(kernlab)
library(RSNNS)
library(caret)
library(ROCR)
library(pROC)
library(ggplot2)
library(dplyr)
library(tidyr)

## set seed
set.seed(1234)

## data
### models
svm_model <- readRDS("/scratch/gtwa/Day/sex_prediction_model/data/celltype_models/models/svm_rad_model_nonneuronal.RDS")
rf_model <- readRDS("/scratch/gtwa/Day/sex_prediction_model/data/celltype_models/models/rf_model_nonneuronal.RDS")
mlp_model <- readRDS("/scratch/gtwa/Day/sex_prediction_model/data/celltype_models/models/mlpwdml_model_nonneuronal.RDS")
lr_model <- readRDS("/scratch/gtwa/Day/sex_prediction_model/data/celltype_models/models/lr_model_nonneuronal.RDS")
chrY_model <- read.csv("/scratch/gtwa/Day/sex_prediction_model/output/celltype_models/Yclass.predictions.nonneuronal.csv", row.names = 1)
Xist_model <- read.csv("/scratch/gtwa/Day/sex_prediction_model/output/celltype_models/Xist.class.predictions.nonneuronal.csv", row.names = 1)

### testing data
testing <- read.csv("/scratch/gtwa/Day/sex_prediction_model/data/celltype_models/test_count_data_nonneuronal.csv", row.names = 1)
#### make character ID bin
testing$Identity_bin_char <- recode(testing$Identity_bin, `0` = "Male", `1` = "Female")
testing$Identity_bin_char <- factor(testing$Identity_bin_char, levels = c("Male", "Female"))
#### make ID bin a factor
testing$Identity_bin <- as.factor(testing$Identity_bin)

# Predictions ####
## logistic regression
lr_prediction <- predict(lr_model,testing[,1:234], type = "response")
## support vector machine
svm_prediction <- predict(svm_model,testing[,1:234], type = "prob")
## random forrest
rf_prediction <- predict(rf_model,testing[,1:234], type = "prob")
## multi-layer perceptron
mlp_prediction <- predict(mlp_model,testing[,1:234], type = "prob")
## chrY expression based prediction
chrY_prediction <- chrY_model$y.counts
## Xist expression based prediction
Xist_prediction <- Xist_model$Xist.counts


# ROC curves ####
## calculate ROC
lr_roc <- roc(testing$Identity_bin,lr_prediction, levels = c(1,0))
svm_roc <- roc(testing$Identity_bin_char, svm_prediction[,2], levels = c("Female", "Male"))
rf_roc <- roc(testing$Identity_bin_char, rf_prediction[,2], levels = c("Female", "Male"))
mlp_roc <- roc(testing$Identity_bin_char, mlp_prediction[,2], levels = c("Female", "Male"))
chrY_roc <- roc(testing$Identity_bin_char, chrY_prediction, levels = c("Female", "Male"))
Xist_roc <- roc(testing$Identity_bin_char, Xist_prediction, levels = c("Female", "Male"), direction = ">")

## plot ROC curves
g <- ggroc(data = list(`Chr Y` = chrY_roc,
                       `Xist` = Xist_roc,
                       `Logistic Regression` = lr_roc,
                       `SVM` = svm_roc,
                       `Random Forest` = rf_roc,
                       `MLP` = mlp_roc)) +
  geom_abline(slope = 1, intercept = 1, linetype = "dashed") +
  labs(color = "Model") +
  theme_bw()
ggsave(filename = "/scratch/gtwa/Day/sex_prediction_model/output/celltype_models/eval/ROC_curves_nonneuronal.pdf",
       plot = g + theme(text = element_text(size = 20)),
       height = 10,
       width = 12,
       device = "pdf")

## save AUC values
auc.table <- data.frame(AUC = c(chrY_roc[["auc"]] %>% as.numeric(),
                                Xist_roc[["auc"]] %>% as.numeric(),
                                lr_roc[["auc"]] %>% as.numeric(),
                                svm_roc[["auc"]] %>% as.numeric(),
                                rf_roc[["auc"]] %>% as.numeric(),
                                mlp_roc[["auc"]] %>% as.numeric()))
rownames(auc.table) <- c("Chr Y", "Xist", "Logistic Regression","SVM", "Random Forest", "MLP")
write.csv(auc.table, "/scratch/gtwa/Day/sex_prediction_model/output/celltype_models/eval/AUCs_table_nonneuronal.csv")

# Confusion Matrices ####
## make precitions using threshold
svm_prediction$class <- ifelse(svm_prediction$Female >0.5, "Female", "Male")
svm_prediction$class <- factor(svm_prediction$class, levels = c("Male", "Female"))

rf_prediction$class <- ifelse(rf_prediction$Female >0.5, "Female", "Male")
rf_prediction$class <- factor(rf_prediction$class, levels = c("Male", "Female"))

mlp_prediction$class <- ifelse(mlp_prediction$Female >0.5, "Female", "Male")
mlp_prediction$class <- factor(mlp_prediction$class, levels = c("Male", "Female"))

lr_prediction <- ifelse(lr_prediction >0.5, 1, 0)
lr_prediction <- as.factor(lr_prediction)

chrY_prediction <- factor(chrY_model$Sex.prediction, levels = c("Male", "Female"))

Xist_prediction <- factor(Xist_model$Sex.prediction, levels = c("Male", "Female"))


## make confusion matrices
svm_confusion <- caret::confusionMatrix(svm_prediction$class,testing$Identity_bin_char)
rf_confusion <- caret::confusionMatrix(rf_prediction$class,testing$Identity_bin_char)
mlp_confusion <- caret::confusionMatrix(mlp_prediction$class,testing$Identity_bin_char)
lr_confusion <- caret::confusionMatrix(lr_prediction, testing$Identity_bin)
chrY_confusion <- caret::confusionMatrix(chrY_prediction, testing$Identity_bin_char)
Xist_confusion <- caret::confusionMatrix(Xist_prediction, testing$Identity_bin_char)

## print confusion matrices to log
print("SVM confusion Matrix: ")
svm_confusion
cat("\n")

print("Random Forest confusion Matrix: ")
rf_confusion
cat("\n")

print("MLP confusion matrix: ")
mlp_confusion
cat("\n")

print("Logistic Regression confusion matrix: ")
lr_confusion
cat("\n")

print("Chr Y confusion matrix: ")
chrY_confusion
cat("\n")

print("Xist confusion matrix: ")
Xist_confusion
cat("\n")

## make accuracy table
acc.table <- data.frame(Accuracy = c(chrY_confusion[["overall"]][["Accuracy"]],
                                     Xist_confusion[["overall"]][["Accuracy"]],
                                     lr_confusion[["overall"]][["Accuracy"]],
                                     svm_confusion[["overall"]][["Accuracy"]],
                                     rf_confusion[["overall"]][["Accuracy"]],
                                     mlp_confusion[["overall"]][["Accuracy"]]))
rownames(acc.table) <- c("Chr Y","Xist", "Logistic Regression", "SVM", "Random Forest", "MLP")
write.csv(acc.table, "/scratch/gtwa/Day/sex_prediction_model/output/celltype_models/eval/Accuracy_table_nonneuronal.csv")

## stacked bar plot of prediction accuracies
### make list of confusion matrices
confusion.list <- list(ChrY = chrY_confusion,
                       Xist = Xist_confusion,
                       `Logistic Regression` = lr_confusion,
                       SVM = svm_confusion,
                       `Random Forest` = rf_confusion,
                       MLP = mlp_confusion)

### extract values to data frame with values for reference_prediction
confusion.df <- data.frame(model = c("Chr Y","Xist", "Logistic Regression", "SVM", "Random Forest", "MLP"),
                           Male_correct = lapply(confusion.list, function(x){x[["table"]][1,1]}) %>% unlist(),
                           Female_incorrect = lapply(confusion.list, function(x){x[["table"]][1,2]}) %>% unlist(),
                           Male_incorrect = lapply(confusion.list, function(x){x[["table"]][2,1]}) %>% unlist(),
                           Female_correct = lapply(confusion.list, function(x){x[["table"]][2,2]}) %>% unlist())

### convert to long format and split reference level and prediction
confusion.df.long <- confusion.df %>%
  pivot_longer(cols = Male_correct:Female_correct,
               names_to = "prediction",
               values_to = "count") %>%
  separate(col = prediction,
           into = c("reference","prediction"),
           sep = "_")

### convert variables to factors
confusion.df.long <- confusion.df.long %>% mutate(reference = factor(reference, levels = c("Male", "Female")),
                                                  prediction = factor(prediction, levels = c("incorrect","correct")),
                                                  model = factor(model, levels = rev(c("Chr Y","Xist", "Logistic Regression", "SVM", "Random Forest", "MLP"))))

### create stacked bar plot of predictions split by sex
confusion.plot <- ggplot(data = confusion.df.long, mapping = aes(x=model, y=count, fill=prediction)) + 
  geom_bar(position="fill", stat="identity") +
  facet_wrap( ~ reference) +
  theme_bw() +
  coord_flip() +
  labs(title = "Confuson plot", x = "Model", y = "Proportion") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

### save
ggsave(plot = confusion.plot,
       file = "/scratch/gtwa/Day/sex_prediction_model/output/celltype_models/eval/confusionplot_nonneuronal.pdf",
       height = 7,
       width = 10,
       device = "pdf")

# sessionInfo
sessionInfo()