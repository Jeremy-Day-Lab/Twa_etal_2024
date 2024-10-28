# Setup ####
## libraries
library(Seurat)
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
svm_model <- readRDS("/path/to/sex_prediction_model/data/models/svm_rad_model.RDS")
rf_model <- readRDS("/path/to/sex_prediction_model/data/models/rf_model.RDS")
mlp_model <- readRDS("/path/to/sex_prediction_model/data/models/mlpwdml_model.RDS")
lr_model <- readRDS("/path/to/sex_prediction_model/data/models/lr_model.RDS")
chrY_model <- read.csv("/path/to/sex_prediction_model/output/Yclass.predictions.csv", row.names = 1)
Xist_model <- read.csv("/path/to/sex_prediction_model/output/Xist.class.predictions.csv", row.names = 1)

### testing data
Rn7_VTA_testing <- readRDS("/path/to/sex_prediction_model/data/Rn7_VTA_testing.RDS")
testing <- read.csv("/path/to/sex_prediction_model/data/test_count_data.csv", row.names = 1)
#### make character ID bin
testing$Identity_bin_char <- recode(testing$Identity_bin, `0` = "Male", `1` = "Female")
testing$Identity_bin_char <- factor(testing$Identity_bin_char, levels = c("Male", "Female"))
#### make ID bin a factor
testing$Identity_bin <- as.factor(testing$Identity_bin)
#### make cell type column
testing$CellType <- Rn7_VTA_testing@meta.data$CellType


# Make predictions ####
## probabilities
### logistic regression
lr_prediction <- predict(lr_model,testing[,1:285], type = "response")
### support vector machine
svm_prediction <- predict(svm_model,testing[,1:285], type = "prob")
### random forrest
rf_prediction <- predict(rf_model,testing[,1:285], type = "prob")
### multi-layer perceptron
mlp_prediction <- predict(mlp_model,testing[,1:285], type = "prob")
### chrY expression based prediction
chrY_prediction <- chrY_model$y.counts
### Xist expression based prediction
Xist_prediction <- Xist_model$Xist.counts

## classify with thresholds
### support vector machine
svm_prediction$class <- ifelse(svm_prediction$Female >0.5, "Female", "Male")
svm_prediction$class <- factor(svm_prediction$class, levels = c("Male", "Female"))
### random forrest
rf_prediction$class <- ifelse(rf_prediction$Female >0.5, "Female", "Male")
rf_prediction$class <- factor(rf_prediction$class, levels = c("Male", "Female"))
### multi-layer perceptron
mlp_prediction$class <- ifelse(mlp_prediction$Female >0.5, "Female", "Male")
mlp_prediction$class <- factor(mlp_prediction$class, levels = c("Male", "Female"))
### logistic regression
lr_prediction <- ifelse(lr_prediction >0.5, 1, 0)
lr_prediction <- as.factor(lr_prediction)
### chrY expression based prediction
chrY_prediction <- factor(chrY_model$Sex.prediction, levels = c("Male", "Female"))
### Xist expression based prediction
Xist_prediction <- factor(Xist_model$Sex.prediction, levels = c("Male", "Female"))


# Create data frame for evaluation ####
eval.df <- data.frame(CellType = testing$CellType,
                      Identity_bin = testing$Identity_bin,
                      Identity_bin_char  = testing$Identity_bin_char,
                      chrY_prediction = chrY_prediction,
                      Xist_prediction = Xist_prediction,
                      lr_prediction = lr_prediction,
                      svm_prediction = svm_prediction$class,
                      rf_prediction = rf_prediction$class,
                      mlp_prediction = mlp_prediction$class)
eval.df <- eval.df %>% mutate(Neuronal_Glial = case_when(CellType %in% c("Glut-Neuron-1", "Glut-Neuron-2", "Glut-Neuron-3", "GABA-Neuron-1", "GABA-Neuron-2", "GABA-Neuron-3", "DA-Neuron") ~ "Neuronal",
                                                         CellType %in% c("Olig-1", "Olig-2", "Olig-3", "Astrocyte", "Polydendrocyte", "Microglia", "OPC-Olig-1", "Mural", "Endothelial") ~ "Glial"))
eval.df.split.celltype <- split(eval.df, f = eval.df$CellType)
eval.df.split.neuronal.glial <- split(eval.df, f = eval.df$Neuronal_Glial)

# Accuracy ####
## Neuronal and Glial cells separately
neuronal.glial.acc <- lapply(eval.df.split.neuronal.glial %>% names(), function(celltype){
  chrY_confusion <- confusionMatrix(eval.df.split.neuronal.glial[[celltype]][["chrY_prediction"]], eval.df.split.neuronal.glial[[celltype]][["Identity_bin_char"]])
  Xist_confusion <- confusionMatrix(eval.df.split.neuronal.glial[[celltype]][["Xist_prediction"]], eval.df.split.neuronal.glial[[celltype]][["Identity_bin_char"]])
  lr_confusion <- confusionMatrix(eval.df.split.neuronal.glial[[celltype]][["lr_prediction"]], eval.df.split.neuronal.glial[[celltype]][["Identity_bin"]])
  svm_confusion <- confusionMatrix(eval.df.split.neuronal.glial[[celltype]][["svm_prediction"]], eval.df.split.neuronal.glial[[celltype]][["Identity_bin_char"]])
  rf_confusion <- confusionMatrix(eval.df.split.neuronal.glial[[celltype]][["rf_prediction"]], eval.df.split.neuronal.glial[[celltype]][["Identity_bin_char"]])
  mlp_confusion <- confusionMatrix(eval.df.split.neuronal.glial[[celltype]][["mlp_prediction"]], eval.df.split.neuronal.glial[[celltype]][["Identity_bin_char"]])
  
  
  ## accuracy table
  acc.table <- data.frame(Accuracy = c(chrY_confusion[["overall"]][["Accuracy"]],
                                       Xist_confusion[["overall"]][["Accuracy"]],
                                       lr_confusion[["overall"]][["Accuracy"]],
                                       svm_confusion[["overall"]][["Accuracy"]],
                                       rf_confusion[["overall"]][["Accuracy"]],
                                       mlp_confusion[["overall"]][["Accuracy"]]))
  rownames(acc.table) <- c("Chr Y", "Xist","Logistic Regression", "SVM", "Random Forest", "MLP")
  
  
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
    labs(title = paste0("Confuson plot: ", celltype), x = "Model", y = "Proportion") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  ### save
  ggsave(plot = confusion.plot,
         file = paste0("/path/to/sex_prediction_model/output/confusionplot.", celltype,".pdf"),
         height = 7,
         width = 10,
         device = "pdf")
  
  # return table
  return(acc.table)
})

names(neuronal.glial.acc) <- names(eval.df.split.neuronal.glial)
#### merge to dataframe 
neuronal.glial.acc.df <- bind_cols(neuronal.glial.acc)
colnames(neuronal.glial.acc.df) <- names(eval.df.split.neuronal.glial)

### Cell Types separately
celltype.acc <- lapply(eval.df.split.celltype %>% names(), function(celltype){
  chrY_confusion <- confusionMatrix(eval.df.split.celltype[[celltype]][["chrY_prediction"]], eval.df.split.celltype[[celltype]][["Identity_bin_char"]])
  Xist_confusion <- confusionMatrix(eval.df.split.celltype[[celltype]][["Xist_prediction"]], eval.df.split.celltype[[celltype]][["Identity_bin_char"]])
  lr_confusion <- confusionMatrix(eval.df.split.celltype[[celltype]][["lr_prediction"]], eval.df.split.celltype[[celltype]][["Identity_bin"]])
  svm_confusion <- confusionMatrix(eval.df.split.celltype[[celltype]][["svm_prediction"]], eval.df.split.celltype[[celltype]][["Identity_bin_char"]])
  rf_confusion <- confusionMatrix(eval.df.split.celltype[[celltype]][["rf_prediction"]], eval.df.split.celltype[[celltype]][["Identity_bin_char"]])
  mlp_confusion <- confusionMatrix(eval.df.split.celltype[[celltype]][["mlp_prediction"]], eval.df.split.celltype[[celltype]][["Identity_bin_char"]])
  
  ## Accuracy table
  acc.table <- data.frame(Accuracy = c(chrY_confusion[["overall"]][["Accuracy"]],
                                       Xist_confusion[["overall"]][["Accuracy"]],
                                       lr_confusion[["overall"]][["Accuracy"]],
                                       svm_confusion[["overall"]][["Accuracy"]],
                                       rf_confusion[["overall"]][["Accuracy"]],
                                       mlp_confusion[["overall"]][["Accuracy"]]))
  rownames(acc.table) <- c("Chr Y", "Xist","Logistic Regression", "SVM", "Random Forest", "MLP")
  
  
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
    labs(title = paste0("Confuson plot: ", celltype), x = "Model", y = "Proportion") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  ### save
  ggsave(plot = confusion.plot,
         file = paste0("/path/to/sex_prediction_model/output/confusionplot.", celltype,".pdf"),
         height = 7,
         width = 10,
         device = "pdf")
  
  
  
  return(acc.table)
})

names(celltype.acc) <- names(eval.df.split.celltype)
#### merge to dataframe
celltype.acc.df <- bind_cols(celltype.acc)
colnames(celltype.acc.df) <- names(eval.df.split.celltype)
#### change column order to be neurons -> glia
celltype.acc.df <- celltype.acc.df %>% select("DA-Neuron", "GABA-Neuron-1", "GABA-Neuron-2", "GABA-Neuron-3", "Glut-Neuron-1", "Glut-Neuron-2", "Glut-Neuron-3", "Olig-1", "Olig-2", "Olig-3", "OPC-Olig-1", "Polydendrocyte", "Astrocyte", "Microglia", "Mural", "Endothelial")


# ROC + AUC ####
## Remake overall table with probabilities
roc.eval.df <- data.frame(CellType = testing$CellType,
                          Identity_bin = testing$Identity_bin,
                          Identity_bin_char  = testing$Identity_bin_char,
                          chrY_prediction = chrY_model$y.counts,
                          Xist_prediction = Xist_model$Xist.counts,
                          lr_prediction = predict(lr_model,testing[,1:285], type = "response"),
                          svm_prediction = svm_prediction[,2],
                          rf_prediction = rf_prediction[,2],
                          mlp_prediction = mlp_prediction[,2])
roc.eval.df <- roc.eval.df %>% mutate(Neuronal_Glial = case_when(CellType %in% c("Glut-Neuron-1", "Glut-Neuron-2", "Glut-Neuron-3", "GABA-Neuron-1", "GABA-Neuron-2", "GABA-Neuron-3", "DA-Neuron") ~ "Neuronal",
                                                                 CellType %in% c("Olig-1", "Olig-2", "Olig-3", "Astrocyte", "Polydendrocyte", "Microglia", "OPC-Olig-1", "Mural", "Endothelial") ~ "Glial"))
roc.eval.df.split.celltype <- split(roc.eval.df, f = roc.eval.df$CellType) # table split by individual cell type
roc.eval.df.split.neuronal.glial <- split(roc.eval.df, f = roc.eval.df$Neuronal_Glial) # table split by neuronal and glial cell types

## Neuronal and glial cell types separately
neuronal.glial.auc <- lapply(roc.eval.df.split.neuronal.glial %>% names(), function(celltype){
  ## calculate ROC
  Xist_roc <- roc(roc.eval.df.split.neuronal.glial[[celltype]][["Identity_bin_char"]],roc.eval.df.split.neuronal.glial[[celltype]][["Xist_prediction"]], levels = c("Female", "Male"), direction = ">")
  chrY_roc <- roc(roc.eval.df.split.neuronal.glial[[celltype]][["Identity_bin_char"]],roc.eval.df.split.neuronal.glial[[celltype]][["chrY_prediction"]], levels = c("Female", "Male"))
  lr_roc <- roc(roc.eval.df.split.neuronal.glial[[celltype]][["Identity_bin"]],roc.eval.df.split.neuronal.glial[[celltype]][["lr_prediction"]], levels = c(1,0))
  svm_roc <- roc(roc.eval.df.split.neuronal.glial[[celltype]][["Identity_bin_char"]],roc.eval.df.split.neuronal.glial[[celltype]][["svm_prediction"]], levels = c("Female", "Male"))
  rf_roc <- roc(roc.eval.df.split.neuronal.glial[[celltype]][["Identity_bin_char"]],roc.eval.df.split.neuronal.glial[[celltype]][["rf_prediction"]], levels = c("Female", "Male"))
  mlp_roc <- roc(roc.eval.df.split.neuronal.glial[[celltype]][["Identity_bin_char"]],roc.eval.df.split.neuronal.glial[[celltype]][["mlp_prediction"]], levels = c("Female", "Male"))
  
  ## ROC curve
  g <- ggroc(data = list(`Chr Y` = chrY_roc,
                         `Xist` = Xist_roc,
                         `Logistic Regression` = lr_roc,
                         `SVM` = svm_roc,
                         `Random Forest` = rf_roc,
                         `MLP` = mlp_roc)) +
    geom_abline(slope = 1, intercept = 1, linetype = "dashed") +
    labs(color = "Model", title = celltype) +
    theme_bw()
  
  ggsave(filename = paste0("/path/to/sex_prediction_model/output/ROC_", celltype,"_curves.pdf"),
         plot = g + theme(text = element_text(size = 20)),
         height = 10,
         width = 12,
         device = "pdf")
  
  ## AUC table
  auc.table <- data.frame(AUC = c(chrY_roc[["auc"]] %>% as.numeric(),
                                  Xist_roc[["auc"]] %>% as.numeric(),
                                  lr_roc[["auc"]] %>% as.numeric(),
                                  svm_roc[["auc"]] %>% as.numeric(),
                                  rf_roc[["auc"]] %>% as.numeric(),
                                  mlp_roc[["auc"]] %>% as.numeric()))
  rownames(auc.table) <- c("Chr Y", "Xist", "Logistic Regression","SVM", "Random Forest", "MLP")
  return(auc.table)
})
names(neuronal.glial.auc) <- names(roc.eval.df.split.neuronal.glial)
### merge to dataframe 
neuronal.glial.auc.df <- bind_cols(neuronal.glial.auc)
colnames(neuronal.glial.auc.df) <- names(roc.eval.df.split.neuronal.glial)



## Cell types separately
celltype.auc <- lapply(roc.eval.df.split.celltype %>% names(), function(celltype){
  ## calculate ROC
  ## calculate ROC
  Xist_roc <- roc(roc.eval.df.split.celltype[[celltype]][["Identity_bin_char"]],roc.eval.df.split.celltype[[celltype]][["Xist_prediction"]], levels = c("Female", "Male"), direction = ">")
  chrY_roc <- roc(roc.eval.df.split.celltype[[celltype]][["Identity_bin_char"]],roc.eval.df.split.celltype[[celltype]][["chrY_prediction"]], levels = c("Female", "Male"))
  lr_roc <- roc(roc.eval.df.split.celltype[[celltype]][["Identity_bin"]],roc.eval.df.split.celltype[[celltype]][["lr_prediction"]], levels = c(1,0))
  svm_roc <- roc(roc.eval.df.split.celltype[[celltype]][["Identity_bin_char"]],roc.eval.df.split.celltype[[celltype]][["svm_prediction"]], levels = c("Female", "Male"))
  rf_roc <- roc(roc.eval.df.split.celltype[[celltype]][["Identity_bin_char"]],roc.eval.df.split.celltype[[celltype]][["rf_prediction"]], levels = c("Female", "Male"))
  mlp_roc <- roc(roc.eval.df.split.celltype[[celltype]][["Identity_bin_char"]],roc.eval.df.split.celltype[[celltype]][["mlp_prediction"]], levels = c("Female", "Male"))
  
  ## ROC curve
  g <- ggroc(data = list(`Chr Y` = chrY_roc,
                         `Xist` = Xist_roc,
                         `Logistic Regression` = lr_roc,
                         `SVM` = svm_roc,
                         `Random Forest` = rf_roc,
                         `MLP` = mlp_roc)) +
    geom_abline(slope = 1, intercept = 1, linetype = "dashed") +
    labs(color = "Model", title = celltype) +
    theme_bw()
  
  ggsave(filename = paste0("/path/to/sex_prediction_model/output/ROC_", celltype,"_curves.pdf"),
         plot = g + theme(text = element_text(size = 20)),
         height = 10,
         width = 12,
         device = "pdf")
  
  ## AUC table
  auc.table <- data.frame(AUC = c(chrY_roc[["auc"]] %>% as.numeric(),
                                  Xist_roc[["auc"]] %>% as.numeric(),
                                  lr_roc[["auc"]] %>% as.numeric(),
                                  svm_roc[["auc"]] %>% as.numeric(),
                                  rf_roc[["auc"]] %>% as.numeric(),
                                  mlp_roc[["auc"]] %>% as.numeric()))
  rownames(auc.table) <- c("Chr Y", "Xist", "Logistic Regression","SVM", "Random Forest", "MLP")
  return(auc.table)
})
names(celltype.auc) <- names(roc.eval.df.split.celltype)
### merge to dataframe
celltype.auc.df <- bind_cols(celltype.auc)
colnames(celltype.auc.df) <- names(roc.eval.df.split.celltype)
### change column order to be neurons -> glia
celltype.auc.df <- celltype.auc.df %>% select("DA-Neuron", "GABA-Neuron-1", "GABA-Neuron-2", "GABA-Neuron-3", "Glut-Neuron-1", "Glut-Neuron-2", "Glut-Neuron-3", "Olig-1", "Olig-2", "Olig-3", "OPC-Olig-1", "Polydendrocyte", "Astrocyte", "Microglia", "Mural", "Endothelial")


# Write out results ####
write.csv(celltype.auc.df, file = "/path/to/sex_prediction_model/output/AUCs_celltype_table.csv")
write.csv(neuronal.glial.auc.df, file = "/path/to/sex_prediction_model/output/AUCs_neuronal_glial_table.csv")
write.csv(celltype.acc.df, file = "/path/to/sex_prediction_model/output/Accuracy_celltype_table.csv")
write.csv(neuronal.glial.acc.df, file = "/path/to/sex_prediction_model/output/Accuracy_neuronal_glial_table.csv")


# Session Info ####
sessionInfo()