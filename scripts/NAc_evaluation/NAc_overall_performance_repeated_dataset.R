# Setup
## libraries
library(Seurat)
library(caret)
library(randomForest)
library(kernlab)
library(RSNNS)
library(ROCR)
library(pROC)
library(ggplot2)
library(dplyr)
library(tidyr)

## set seed
set.seed(1234)

## functions
### Chr Y classification
y.count <- function(srat.object){
  # load gene expression count data
  srat.counts <- as.data.frame(t(as.matrix(GetAssayData(object = srat.object, slot = "data", assay = "RNA"))))
  
  # load gene annotations
  Rn7.gtf <- rtracklayer::import("/scratch/gtwa/Day/sex_prediction_model/data/Rattus_norvegicus.mRatBN7.2.105.gtf")
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
    rowSums() %>% data.frame()
  ### set column name
  colnames(y.counts) <- c("y.counts")
  
  return(y.counts)
}
y.count.classify <- function(srat.object, threshold){
  srat.object$y.counts <- y.count(srat.object)
  classifications <- case_when(srat.object$y.counts > threshold ~ "Male",
                               TRUE ~ "Female")
  return(classifications)
}

### Xist classification
Xist.count <- function(srat.object){
  # load gene expression count data
  Xist.counts <- as.data.frame(t(as.matrix(GetAssayData(object = srat.object, slot = "data", assay = "RNA")))) %>%
    pull("ENSRNOG00000065796") # pull out Xist (named as ENSRNOG00000065796)
  
  return(Xist.counts)
}
Xist.count.classify <- function(srat.object, threshold){
  srat.object$Xist.counts <- Xist.count(srat.object)
  classifications <- case_when(srat.object$Xist.counts > threshold ~ "Female",
                               TRUE ~ "Male")
  return(classifications)
}

## data
### NAc object
testing.NAc <- readRDS("/data/project/daylab/2019-JD-0040/MCN_Code/Objects/NAc_Combo_Integrated.RDS")
#### subset for repeated data set
testing.NAc@active.ident <- factor(testing.NAc$Dataset)
testing.NAc <- subset(testing.NAc, idents = "Repeated")

### model genes
confirmed.genes <- read.csv("/scratch/gtwa/Day/sex_prediction_model/data/Boruta_final_decisions.csv") %>% 
  filter(finalDecision == "Confirmed") %>%
  pull(variable)
### models
svm_model <- readRDS("/scratch/gtwa/Day/sex_prediction_model/data/models/svm_rad_model.RDS")
rf_model <- readRDS("/scratch/gtwa/Day/sex_prediction_model/data/models/rf_model.RDS")
mlp_model <- readRDS("/scratch/gtwa/Day/sex_prediction_model/data/models/mlpwdml_model.RDS")
lr_model <- readRDS("/scratch/gtwa/Day/sex_prediction_model/data/models/lr_model.RDS")

# NAc object preparation ####
## extract count matrices
testing.NAc.df <- as.data.frame(t(as.matrix(GetAssayData(object = testing.NAc, slot = "data",assay = "RNA")))) %>% 
  select(any_of(confirmed.genes))
### make sure cells are in the same order
all(names(Idents(testing.NAc)) == rownames(testing.NAc.df)) #TRUE 
### replace gene names that have "-" with "." to match gene name form the VTA models expect
colnames(testing.NAc.df) <- gsub("-", ".", colnames(testing.NAc.df))

## add sex identity columns
### character column
testing.NAc.df$Identity_bin_char <- factor(testing.NAc$Sex, levels = c("Male", "Female"))
### binary columns 
testing.NAc.df$Identity_bin <- ifelse(testing.NAc.df$Identity_bin_char  == "Female",
                                      1, #Females are 1
                                      0) #Males are 0
testing.NAc.df$Identity_bin <- as.factor(testing.NAc.df$Identity_bin)

# Make predictions ####
## logistic regression
lr_prediction <- predict(lr_model,testing.NAc.df[,1:285], type = "response")
## support vector machine
svm_prediction <- predict(svm_model,testing.NAc.df[,1:285], type = "prob")
## random forest
rf_prediction <- predict(rf_model,testing.NAc.df[,1:285], type = "prob")
## multi-layer perceptron
mlp_prediction <- predict(mlp_model,testing.NAc.df[,1:285], type = "prob")
## chrY expression based prediction
chrY_prediction <- data.frame(cell = names(Idents(testing.NAc)),
                              sex = testing.NAc$Sex,
                              counts = y.count(testing.NAc),
                              prediction = y.count.classify(testing.NAc,0))
chrY_prediction$prediction <- factor(chrY_prediction$prediction, levels = c("Male", "Female"))
## Xist expression based prediction
Xist_prediction <- data.frame(cell = names(Idents(testing.NAc)),
                              sex = testing.NAc$Sex,
                              counts = Xist.count(testing.NAc),
                              prediction = Xist.count.classify(testing.NAc,0))
Xist_prediction$prediction <- factor(Xist_prediction$prediction, levels = c("Male", "Female"))


# ROC + AUC ####
lr_roc <- roc(testing.NAc.df$Identity_bin,lr_prediction, levels = c(1,0))
svm_roc <- roc(testing.NAc.df$Identity_bin_char, svm_prediction[,2], levels = c("Female", "Male"))
rf_roc <- roc(testing.NAc.df$Identity_bin_char, rf_prediction[,2], levels = c("Female", "Male"))
mlp_roc <- roc(testing.NAc.df$Identity_bin_char, mlp_prediction[,2], levels = c("Female", "Male"))
chrY_roc <- roc(testing.NAc.df$Identity_bin_char, chrY_prediction$y.counts, levels = c("Female", "Male"))
Xist_roc <- roc(testing.NAc.df$Identity_bin_char, Xist_prediction$counts, levels = c("Female", "Male"), direction = ">")

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
ggsave(filename = "/scratch/gtwa/Day/sex_prediction_model/output/NAc_evaluation/ROC_curves_repeated.pdf",
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
write.csv(auc.table, "/scratch/gtwa/Day/sex_prediction_model/output/NAc_evaluation/AUCs_table_repeated.csv")

# Confusion matrices ####
## make precitions using threshold
svm_prediction$class <- ifelse(svm_prediction$Female >0.5, "Female", "Male")
svm_prediction$class <- factor(svm_prediction$class, levels = c("Male", "Female"))

rf_prediction$class <- ifelse(rf_prediction$Female >0.5, "Female", "Male")
rf_prediction$class <- factor(rf_prediction$class, levels = c("Male", "Female"))

mlp_prediction$class <- ifelse(mlp_prediction$Female >0.5, "Female", "Male")
mlp_prediction$class <- factor(mlp_prediction$class, levels = c("Male", "Female"))

lr_prediction <- ifelse(lr_prediction >0.5, 1, 0)
lr_prediction <- as.factor(lr_prediction)

## make confusion matrices
svm_confusion <- caret::confusionMatrix(svm_prediction$class, testing.NAc.df$Identity_bin_char)
rf_confusion <- caret::confusionMatrix(rf_prediction$class, testing.NAc.df$Identity_bin_char)
mlp_confusion <- caret::confusionMatrix(mlp_prediction$class, testing.NAc.df$Identity_bin_char)
lr_confusion <- caret::confusionMatrix(lr_prediction, testing.NAc.df$Identity_bin)
chrY_confusion <- caret::confusionMatrix(chrY_prediction$prediction, testing.NAc.df$Identity_bin_char)
Xist_confusion <- caret::confusionMatrix(Xist_prediction$prediction, testing.NAc.df$Identity_bin_char)

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
write.csv(acc.table, "/scratch/gtwa/Day/sex_prediction_model/output/NAc_evaluation/Accuracy_table_repeated.csv")

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
       file = "/scratch/gtwa/Day/sex_prediction_model/output/NAc_evaluation/confusionplot_repeated.pdf",
       height = 7,
       width = 10,
       device = "pdf")

# sessionInfo ####
sessionInfo()