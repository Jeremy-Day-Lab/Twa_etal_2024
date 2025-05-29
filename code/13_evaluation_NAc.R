# Setup ####
# libraries
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

# set seet
set.seed(1234)

# Data
## gene chromosome data
Rn7.gtf <- rtracklayer::import("/scratch/gtwa/Day/sex_prediction_model/Rn7_Xist_add/raw_data/Rattus_norvegicus.mRatBN7.2.105.Xist.gtf")
Rn7.gtf <- as.data.frame(Rn7.gtf)

# Process NAc data ####
# Load NAc data
## NAc w/o Xist
NAc_noXist <- readRDS("/scratch/gtwa/Day/sex_prediction_model/Rn7_Xist_add/raw_data/NAc_Combo_Integrated.RDS")
## Acute data w/ Xist
acute <- list(
  Fem_Sal  = Read10X(data.dir = "/data/project/daylab/2019-JD-0037/Ensembl_Rn7_Xist/FemSal_output_Rn7_Xist/outs/filtered_feature_bc_matrix/"),
  Fem_Coc  = Read10X(data.dir = "/data/project/daylab/2019-JD-0037/Ensembl_Rn7_Xist/FemCoc_output_Rn7_Xist/outs/filtered_feature_bc_matrix/"),
  Male_Sal = Read10X(data.dir = "/data/project/daylab/2019-JD-0037/Ensembl_Rn7_Xist/MaleSal_output_Rn7_Xist/outs/filtered_feature_bc_matrix/"),
  Male_Coc = Read10X(data.dir = "/data/project/daylab/2019-JD-0037/Ensembl_Rn7_Xist/MaleCoc_output_Rn7_Xist/outs/filtered_feature_bc_matrix/")
)
## Repeated data w/ Xist
repeated <- list(
  Fem_Sal  = Read10X(data.dir = "/data/project/daylab/2019-JD-0040/Ensembl_Rn7_Xist/1_1_output_Rn7_Xist/outs/filtered_feature_bc_matrix/"),
  Fem_Coc  = Read10X(data.dir = "/data/project/daylab/2019-JD-0040/Ensembl_Rn7_Xist/2_1_output_Rn7_Xist/outs/filtered_feature_bc_matrix/"),
  Male_Sal = Read10X(data.dir = "/data/project/daylab/2019-JD-0040/Ensembl_Rn7_Xist/3_1_output_Rn7_Xist/outs/filtered_feature_bc_matrix/"),
  Male_Coc = Read10X(data.dir = "/data/project/daylab/2019-JD-0040/Ensembl_Rn7_Xist/4_1_output_Rn7_Xist/outs/filtered_feature_bc_matrix/")
)

# Create Seurat objects
acute <- lapply(acute, function(x) {
  CreateSeuratObject(counts = x, project = "Acute", min.cells = 0, min.features = 0)
})
repeated <- lapply(repeated, function(x) {
  CreateSeuratObject(counts = x, project = "Repeated", min.cells = 0, min.features = 0)
})

# Add metadata
acute_names <- names(acute)
acute <- lapply(names(acute), function(x) {
  acute[[x]]$Dataset <- "Acute"
  acute[[x]]$Sex_Stim <- x
  return(acute[[x]])
})
names(acute) <- acute_names

repeated_names <- names(repeated)
repeated <- lapply(names(repeated), function(x) {
  repeated[[x]]$Dataset <- "Repeated"
  repeated[[x]]$Sex_Stim <- x
  return(repeated[[x]])
})
names(repeated) <- repeated_names

# determine the appropriate sample tags from the !Xist object so cells can be matched
extract_tags <- function(metadata) {
  # Extract the tag from the cell names
  metadata$Tag <- sub(".*-1", "", rownames(metadata))
  
  # Group by Dataset and Sex_Stim and get the unique tag for each combination
  result <- metadata %>%
    dplyr::group_by(Dataset, Sex_Stim) %>%
    dplyr::summarize(Tag = unique(Tag), .groups = "drop")
  
  return(result)
}

sample_tags <- extract_tags(NAc_noXist@meta.data)
# Dataset  Sex_Stim Tag  
# 1 Acute    Fem_Coc  _2_1 
# 2 Acute    Fem_Sal  _1_1 
# 3 Acute    Male_Coc _4_1 
# 4 Acute    Male_Sal _3_1 
# 5 Repeated Fem_Coc  _2_2 
# 6 Repeated Fem_Sal  _1_2 
# 7 Repeated Male_Coc _4_2 
# 8 Repeated Male_Sal _3_2 

# Append tags to acute and repeated data cell names
acute <- lapply(acute, function(x) {
  # Get the corresponding tag from the sample_tags data frame
  tag <- sample_tags$Tag[sample_tags$Dataset == x$Dataset[1] & sample_tags$Sex_Stim == x$Sex_Stim[1]]
  
  # Append the tag to the cell names
  x <- RenameCells(x, new.names = paste0(colnames(x), tag))
  return(x)
})
repeated <- lapply(repeated, function(x) {
  # Get the corresponding tag from the sample_tags data frame
  tag <- sample_tags$Tag[sample_tags$Dataset == x$Dataset[1] & sample_tags$Sex_Stim == x$Sex_Stim[1]]
  
  # Append the tag to the cell names
  x <- RenameCells(x, new.names = paste0(colnames(x), tag))
  return(x)
})

# Verify that the tags were added correctly, and match a subset of the orignial data
acute_cells <- lapply(acute, function(x) {
  matching_cells <- NAc_noXist@meta.data %>%
    filter(Dataset == x$Dataset[1] & Sex_Stim == x$Sex_Stim[1]) %>%
    rownames()
  cells <- Cells(x)[Cells(x) %in% matching_cells]
  return(cells)
})
repeated_cells <- lapply(repeated, function(x) {
  matching_cells <- NAc_noXist@meta.data %>%
    filter(Dataset == x$Dataset[1] & Sex_Stim == x$Sex_Stim[1]) %>%
    rownames()
  cells <- Cells(x)[Cells(x) %in% matching_cells]
  return(cells)
})
## create a list of all cells that are shared between Rn7_Xist and the Xist added data
shared_cells <- c(unlist(acute_cells), unlist(repeated_cells))

# Extract raw Xist counts from the acute and repeated data
Xist_acute <- lapply(names(acute), FUN = function(x,y) {
  xist_counts <- LayerData(acute[[x]], "counts", cells = acute_cells[[x]], features = "Xist")
  return(xist_counts)
})
names(Xist_acute) <- acute_names

Xist_repeated <- lapply(names(repeated), function(x) {
  xist_counts <- LayerData(repeated[[x]], "counts", cells = repeated_cells[[x]], features = "Xist")
  return(xist_counts)
})
names(Xist_repeated) <- repeated_names

## Combine acute and repeated data
Xist_counts <- do.call(cbind, c(Xist_acute, Xist_repeated))

# Extract raw counts from NAc_noXist data
NAc_noXist_counts <- LayerData(NAc_noXist, "counts", cells = shared_cells)
## Merge Xist counts with NAc_noXist counts
NAc_Xist_counts <- rbind(NAc_noXist_counts, Xist_counts)

# Create a new Seurat object with the merged counts
NAc_Xist <- CreateSeuratObject(counts = NAc_Xist_counts, project = "NAc_Xist")
## Add metadata
NAc_Xist$Dataset <- NAc_noXist$Dataset[rownames(NAc_Xist@meta.data)]
NAc_Xist$Sex <- NAc_noXist$Sex[rownames(NAc_Xist@meta.data)]
NAc_Xist$Stim <- NAc_noXist$Stim[rownames(NAc_Xist@meta.data)]
NAc_Xist$Sex_Stim <- NAc_noXist$Sex_Stim[rownames(NAc_Xist@meta.data)]
NAc_Xist$CellType <- NAc_noXist$CellType[rownames(NAc_Xist@meta.data)]
## Normalise data
NAc_Xist <- NormalizeData(NAc_Xist, normalization.method = "LogNormalize", scale.factor = 10000)
## Scale data
NAc_Xist <- ScaleData(NAc_Xist, verbose = FALSE)
## Save object
saveRDS(NAc_Xist, file = "/scratch/gtwa/Day/sex_prediction_model/Rn7_Xist_add/processed_data/13_evaluation_NAc/NAc_Xist.RDS")


# Evaluate Model Performance ####
# Trained models
lr_model <- readRDS("/scratch/gtwa/Day/sex_prediction_model/Rn7_Xist_add/processed_data/5_training_logistic_regression/lr_model.RDS")
rf_model <- readRDS("/scratch/gtwa/Day/sex_prediction_model/Rn7_Xist_add/processed_data/6_training_random_forest/rf_model.RDS")
svm_model <- readRDS("/scratch/gtwa/Day/sex_prediction_model/Rn7_Xist_add/processed_data/7_training_svm/svm_rad_model.RDS")
mlp_model <- readRDS("/scratch/gtwa/Day/sex_prediction_model/Rn7_Xist_add/processed_data/8_training_mlp/mlpwdml_model.RDS")

# Baseline models
## Chr Y
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

# Create count data subset for models
confirmed.genes <- read.csv("/scratch/gtwa/Day/sex_prediction_model/Rn7_Xist_add/processed_data/4_feature_selection_Boruta/Boruta_final_decisions.csv") %>% 
  filter(finalDecision == "Confirmed") %>%
  pull(variable)

NAc_df <- LayerData(NAc_Xist, "data", features = confirmed.genes) %>%
  as.matrix() %>%
  t() %>%
  as.data.frame() %>%
  mutate(Identity_bin_char = factor(NAc_Xist$Sex, levels = c("Male", "Female")),
         Identity_bin = Identity_bin_char %>% as.numeric() - 1)

## replace "-" in gene names with "."
colnames(NAc_df) <- gsub("-", ".", colnames(NAc_df))

## Predictions ####
# logistic regression
lr_prediction <- predict(lr_model,NAc_df[,1:(ncol(NAc_df)-2)], type = "response")
# support vector machine
svm_prediction <- predict(svm_model,NAc_df[,1:(ncol(NAc_df)-2)], type = "prob")
# random forrest
rf_prediction <- predict(rf_model,NAc_df[,1:(ncol(NAc_df)-2)], type = "prob")
# multi-layer perceptron
mlp_prediction <- predict(mlp_model,NAc_df[,1:(ncol(NAc_df)-2)], type = "prob")
# chrY expression based prediction
chrY_prediction <- y.count(NAc_Xist)
# Xist expression based prediction
Xist_prediction <- Xist.count(NAc_Xist)


## ROC curves ####
# calculate ROC
lr_roc <- roc(NAc_df$Identity_bin,lr_prediction, levels = c(1,0))
svm_roc <- roc(NAc_df$Identity_bin_char, svm_prediction[,2], levels = c("Female", "Male"))
rf_roc <- roc(NAc_df$Identity_bin_char, rf_prediction[,2], levels = c("Female", "Male"))
mlp_roc <- roc(NAc_df$Identity_bin_char, mlp_prediction[,2], levels = c("Female", "Male"))
chrY_roc <- roc(NAc_df$Identity_bin_char, chrY_prediction, levels = c("Female", "Male"))
Xist_roc <- roc(NAc_df$Identity_bin_char, Xist_prediction, levels = c("Female", "Male"), direction = ">")

# plot ROC curves
g <- ggroc(data = list(`Chr Y` = chrY_roc,
                       `Xist` = Xist_roc,
                       `Logistic Regression` = lr_roc,
                       `SVM` = svm_roc,
                       `Random Forest` = rf_roc,
                       `MLP` = mlp_roc), legacy.axes = TRUE) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(color = "Model") +
  theme_bw()
ggsave(filename = "/data/scratch/gtwa/Day/sex_prediction_model/Rn7_Xist_add/plots/13_evaluation_NAc/ROC_curves.pdf",
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
rownames(auc.table) <- c("Chr Y",
                         "Xist",
                         "Logistic Regression",
                         "SVM",
                         "Random Forest",
                         "MLP")
write.csv(auc.table, "/data/scratch/gtwa/Day/sex_prediction_model/Rn7_Xist_add/processed_data/13_evaluation_NAc/AUCs_table.csv")

## Confusion Matrices ####
# make precitions using threshold
svm_prediction$class <- ifelse(svm_prediction$Female >0.5, "Female", "Male")
svm_prediction$class <- factor(svm_prediction$class, levels = c("Male", "Female"))

rf_prediction$class <- ifelse(rf_prediction$Female >0.5, "Female", "Male")
rf_prediction$class <- factor(rf_prediction$class, levels = c("Male", "Female"))

mlp_prediction$class <- ifelse(mlp_prediction$Female >0.5, "Female", "Male")
mlp_prediction$class <- factor(mlp_prediction$class, levels = c("Male", "Female"))

lr_prediction <- ifelse(lr_prediction >0.5, 1, 0)
lr_prediction <- as.factor(lr_prediction)

chrY_class <- ifelse(chrY_prediction > 0.5, "Male", "Female")
chrY_class <- factor(chrY_class, levels = c("Male", "Female"))

Xist_class <- ifelse(Xist_prediction < 0.5, "Male", "Female")
Xist_class <- factor(Xist_class, levels = c("Male", "Female"))


# make confusion matrices
svm_confusion <- caret::confusionMatrix(svm_prediction$class,NAc_df$Identity_bin_char)
rf_confusion <- caret::confusionMatrix(rf_prediction$class,NAc_df$Identity_bin_char)
mlp_confusion <- caret::confusionMatrix(mlp_prediction$class,NAc_df$Identity_bin_char)
lr_confusion <- caret::confusionMatrix(lr_prediction, as.factor(NAc_df$Identity_bin))
chrY_confusion <- caret::confusionMatrix(chrY_class, NAc_df$Identity_bin_char)
Xist_confusion <- caret::confusionMatrix(Xist_class, NAc_df$Identity_bin_char)

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
rownames(acc.table) <- c("Chr Y",
                         "Xist",
                         "Logistic Regression",
                         "SVM",
                         "Random Forest",
                         "MLP")
write.csv(acc.table, "/data/scratch/gtwa/Day/sex_prediction_model/Rn7_Xist_add/processed_data/13_evaluation_NAc/Accuracy_table.csv")

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
  scale_fill_manual(values = c(correct = "#00A651", incorrect = "#000000")) +
  geom_bar(position="fill", stat="identity") +
  facet_wrap( ~ reference) +
  theme_bw() +
  coord_flip() +
  labs(title = "Confuson plot", x = "Model", y = "Proportion") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

### save
ggsave(plot = confusion.plot,
       file = "/data/scratch/gtwa/Day/sex_prediction_model/Rn7_Xist_add/plots/13_evaluation_NAc/confusionplot.pdf",
       height = 7,
       width = 10,
       device = "pdf")

# sessionInfo
sessionInfo()
