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

### models
lr_model <- readRDS("/scratch/gtwa/Day/sex_prediction_model/Rn7_Xist_add/processed_data/5_training_logistic_regression/lr_model.RDS")
rf_model <- readRDS("/scratch/gtwa/Day/sex_prediction_model/Rn7_Xist_add/processed_data/6_training_random_forest/rf_model.RDS")
svm_model <- readRDS("/scratch/gtwa/Day/sex_prediction_model/Rn7_Xist_add/processed_data/7_training_svm/svm_rad_model.RDS")
mlp_model <- readRDS("/scratch/gtwa/Day/sex_prediction_model/Rn7_Xist_add/processed_data/8_training_mlp/mlpwdml_model.RDS")

### testing data
#### testing object
testing_obj <- readRDS("/scratch/gtwa/Day/sex_prediction_model/Rn7_Xist_add/processed_data/2_create_splits/Rn7_VTA_testing.RDS")
#### count matrix
testing <- read.csv("/scratch/gtwa/Day/sex_prediction_model/Rn7_Xist_add/processed_data/4_feature_selection_Boruta/test_count_data.csv", row.names = 1)
#### make character ID bin
testing$Identity_bin_char <- recode(testing$Identity_bin, `0` = "Male", `1` = "Female")
testing$Identity_bin_char <- factor(testing$Identity_bin_char, levels = c("Male", "Female"))
#### make ID bin a factor
testing$Identity_bin <- as.factor(testing$Identity_bin)
#### make cell type column
testing$CellType <- testing_obj@meta.data$CellType


# Baseline Classifiers ####
## chr Y
y.count <- function(srat.object){
  # load gene expression count data
  srat.counts <- as.data.frame(t(as.matrix(GetAssayData(object = srat.object, slot = "data", assay = "RNA"))))
  
  # load gene annotations
  Rn7.gtf <- rtracklayer::import("/scratch/gtwa/Day/sex_prediction_model/Rn7_Xist_add/raw_data/Rattus_norvegicus.mRatBN7.2.105.Xist.gtf")
  Rn7.gtf <- as.data.frame(Rn7.gtf)
  
  # select Y chromosome genes
  y.genes <- Rn7.gtf %>% #from Rn7 gene annotations
    filter(seqnames == "Y") %>% #filter for chrY only
    select(gene_id, gene_name) %>% #select gene ids and names
    mutate(gene_name = coalesce(gene_name, gene_id)) %>% #replace NA name values with gene ID
    pull(gene_name) %>% #pull gene name values
    unique() #only unique values, returns 27 genes
  
  # create cumulative expresion of Y genes
  y.counts <- srat.counts %>%
    select(any_of(y.genes)) %>%
    rowSums()
  
  return(y.counts)
}
## Xist
Xist.count <- function(srat.object){
  # load gene expression count data
  Xist.counts <- as.data.frame(t(as.matrix(GetAssayData(object = srat.object, slot = "data", assay = "RNA")))) %>%
    pull("Xist")
  
  return(Xist.counts)
}

# Predictions ####
# class probabilities
## logistic regression
lr_prediction <- predict(lr_model,testing[,1:(ncol(testing)-3)], type = "response")
## support vector machine
svm_prediction <- predict(svm_model,testing[,1:(ncol(testing)-3)], type = "prob")
## random forest
rf_prediction <- predict(rf_model,testing[,1:(ncol(testing)-3)], type = "prob")
## multi-layer perceptron
mlp_prediction <- predict(mlp_model,testing[,1:(ncol(testing)-3)], type = "prob")
## chrY expression based prediction
chrY_prediction <- y.count(testing_obj)
## Xist expression based prediction
Xist_prediction <- Xist.count(testing_obj)

# binary classification
## logistic regression
lr_class <- ifelse(lr_prediction >0.5, 1, 0)
lr_class <- as.factor(lr_class)
## support vector machine
svm_prediction$class <- ifelse(svm_prediction$Female >0.5, "Female", "Male")
svm_prediction$class <- factor(svm_prediction$class, levels = c("Male", "Female"))
## random forest
rf_prediction$class <- ifelse(rf_prediction$Female >0.5, "Female", "Male")
rf_prediction$class <- factor(rf_prediction$class, levels = c("Male", "Female"))
## multi-layer perceptron
mlp_prediction$class <- ifelse(mlp_prediction$Female >0.5, "Female", "Male")
mlp_prediction$class <- factor(mlp_prediction$class, levels = c("Male", "Female"))
## chrY expression based prediction
chrY_class <- ifelse(chrY_prediction > 0.5, "Male", "Female")
chrY_class <- factor(chrY_class, levels = c("Male", "Female"))
## Xist expression based prediction
Xist_class <- ifelse(Xist_prediction < 0.5, "Male", "Female")
Xist_class <- factor(Xist_class, levels = c("Male", "Female"))

# Create data frame for evaluation ####
eval.df <- data.frame(CellType = testing$CellType,
                      Identity_bin = testing$Identity_bin,
                      Identity_bin_char  = testing$Identity_bin_char,
                      chrY_prediction = chrY_class,
                      Xist_prediction = Xist_class,
                      lr_prediction = lr_class,
                      svm_prediction = svm_prediction$class,
                      rf_prediction = rf_prediction$class,
                      mlp_prediction = mlp_prediction$class)
eval.df <- eval.df %>% mutate(Neuronal_Glial = case_when(CellType %in% c("Glut-Neuron-1", "Glut-Neuron-2", "Glut-Neuron-3", "GABA-Neuron-1", "GABA-Neuron-2", "GABA-Neuron-3", "DA-Neuron") ~ "Neuronal",
                                                         CellType %in% c("Olig-1", "Olig-2", "Olig-3", "Astrocyte", "Polydendrocyte", "Microglia", "OPC-Olig-1", "Mural", "Endothelial") ~ "Glial"))
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
         file = paste0("/scratch/gtwa/Day/sex_prediction_model/Rn7_Xist_add/plots/10_evaluation_celltype_specific/confusionplot.", celltype,".pdf"),
         height = 7,
         width = 10,
         device = "pdf")
  
  # return table
  return(acc.table)
})

names(neuronal.glial.acc) <- names(eval.df.split.neuronal.glial)
# merge to dataframe 
neuronal.glial.acc.df <- bind_cols(neuronal.glial.acc)
colnames(neuronal.glial.acc.df) <- names(eval.df.split.neuronal.glial)
## write out
write.csv(neuronal.glial.acc.df, file = "/scratch/gtwa/Day/sex_prediction_model/Rn7_Xist_add/processed_data/10_evaluation_celltype_specific/Accuracy_neuronal_glial_table.csv")

# ROC + AUC ####
roc.eval.df <- data.frame(CellType = testing$CellType,
                          Identity_bin = testing$Identity_bin,
                          Identity_bin_char  = testing$Identity_bin_char,
                          chrY_prediction = chrY_prediction,
                          Xist_prediction = Xist_prediction,
                          lr_prediction = lr_prediction,
                          svm_prediction = svm_prediction[,2],
                          rf_prediction = rf_prediction[,2],
                          mlp_prediction = mlp_prediction[,2])
roc.eval.df <- roc.eval.df %>% mutate(Neuronal_Glial = case_when(CellType %in% c("Glut-Neuron-1", "Glut-Neuron-2", "Glut-Neuron-3", "GABA-Neuron-1", "GABA-Neuron-2", "GABA-Neuron-3", "DA-Neuron") ~ "Neuronal",
                                                                 CellType %in% c("Olig-1", "Olig-2", "Olig-3", "Astrocyte", "Polydendrocyte", "Microglia", "OPC-Olig-1", "Mural", "Endothelial") ~ "Glial"))
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
  
  ggsave(filename = paste0("/scratch/gtwa/Day/sex_prediction_model/Rn7_Xist_add/plots/10_evaluation_celltype_specific/ROC_", celltype,"_curves.pdf"),
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
# Write out results
write.csv(neuronal.glial.auc.df, file = "/scratch/gtwa/Day/sex_prediction_model/Rn7_Xist_add/processed_data/10_evaluation_celltype_specific/AUCs_neuronal_glial_table.csv")

# sesionInfo ####
sessionInfo()