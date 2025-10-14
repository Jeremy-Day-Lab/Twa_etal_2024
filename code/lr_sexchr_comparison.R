# Setup ####
## libraries
library(Seurat)
library(caret)
library(ROCR)
library(pROC)
library(ggplot2)
library(dplyr)
library(tidyr)

## set seed
set.seed(1234)

## data
### gene chromosome data
Rn7.gtf <- rtracklayer::import("/scratch/gtwa/Day/sex_prediction_model/Rn7_Xist_add/raw_data/Rattus_norvegicus.mRatBN7.2.105.Xist.gtf")
Rn7.gtf <- as.data.frame(Rn7.gtf)

### training data
#### training object
training_obj <- readRDS("/scratch/gtwa/Day/sex_prediction_model/Rn7_Xist_add/processed_data/2_create_splits/Rn7_VTA_training.RDS")

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
testing <- testing %>% mutate(CellType = case_when(CellType %in% c("Glut-Neuron-1", "Glut-Neuron-2", "Glut-Neuron-3", "GABA-Neuron-1", "GABA-Neuron-2", "GABA-Neuron-3", "DA-Neuron") ~ "Neuronal",
                                                   CellType %in% c("Olig-1", "Olig-2", "Olig-3", "Astrocyte", "Polydendrocyte", "Microglia", "OPC-Olig-1", "Mural", "Endothelial") ~ "NonNeuronal"))

### boruta Feature Selection
boruta_decisions <- read.csv("/scratch/gtwa/Day/sex_prediction_model/Rn7_Xist_add/processed_data/4_feature_selection_Boruta/Boruta_final_decisions.csv")
boruta_features <- boruta_decisions$variable[boruta_decisions$finalDecision == "Confirmed"]


# Select model features ####
# Chr Y Genes
chrY_genes <- function(srat.object){
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
    unique() #only unique values
  
  # select gene names
  genes <- srat.counts %>%
    select(any_of(y.genes)) %>%
    colnames()
  
  return(genes)
}
y_genes <- chrY_genes(training_obj)
# Chr X Genes
chrX_genes <- function(srat.object){
  # load gene expression count data
  srat.counts <- as.data.frame(t(as.matrix(GetAssayData(object = srat.object, slot = "data", assay = "RNA"))))
  
  # load gene annotations
  Rn7.gtf <- rtracklayer::import("/scratch/gtwa/Day/sex_prediction_model/Rn7_Xist_add/raw_data/Rattus_norvegicus.mRatBN7.2.105.Xist.gtf")
  Rn7.gtf <- as.data.frame(Rn7.gtf)
  
  # select X chromosome genes
  x.genes <- Rn7.gtf %>% #from Rn7 gene annotations
    filter(seqnames == "X") %>% #filter for chrX only
    select(gene_id, gene_name) %>% #select gene ids and names
    mutate(gene_name = coalesce(gene_name, gene_id)) %>% #replace NA name values with gene ID
    pull(gene_name) %>% #pull gene name values
    unique() #only unique values
  
  # select gene names
  genes <- srat.counts %>%
    select(any_of(x.genes)) %>%
    colnames()
  
  return(genes)
}
x_genes <- chrX_genes(training_obj)
# Sex Chromosome genes
sexchr_features <- boruta_features[boruta_features %in% c(y_genes, x_genes)]
print("Sex Chromosome genes included in Boruta: ")
print(sexchr_features)
cat("\n")
# boruta features - (Chr Y genes + Xist)
non_sexchr_features <- boruta_features[!boruta_features %in% c(y_genes, x_genes)]


# Prepare training data ####
# Sex chromosome features only
training_sexchr <- as.data.frame(t(as.matrix(GetAssayData(object = training_obj,slot = "data",assay = "RNA")))) %>% 
  select(any_of(sexchr_features))
training_sexchr$Identity_bin <- ifelse(training_obj$Sex == "Female", 1, 0)

# Non Sex chromosome features only
training_non_sexchr <- as.data.frame(t(as.matrix(GetAssayData(object = training_obj,slot = "data",assay = "RNA")))) %>% 
  select(any_of(non_sexchr_features))
training_non_sexchr$Identity_bin <- ifelse(training_obj$Sex == "Female", 1, 0)
## replace "-" in gene names with "."
colnames(training_non_sexchr) <- gsub("-", ".", colnames(training_non_sexchr))

# LR Models ####
## All features LR model ####
lr_model_all <- readRDS("/scratch/gtwa/Day/sex_prediction_model/Rn7_Xist_add/processed_data/5_training_logistic_regression/lr_model.RDS")

## New models ####
# sex chromosome features only
lr_model_sexchr <- glm(Identity_bin ~ .,
                       data = training_sexchr,
                       family = "binomial")
# model summary
summary(lr_model_sexchr)
# save model
saveRDS(lr_model_sexchr, file = "/scratch/gtwa/Day/sex_prediction_model/Rn7_Xist_add/processed_data/lr_sexchr_comparison/sexchr_only_lr.RDS")

# non sex chromosome features only
lr_model_nonsexchr <- glm(Identity_bin ~ .,
                          data = training_non_sexchr,
                          family = "binomial")
# model summary
summary(lr_model_nonsexchr)
# save model
saveRDS(lr_model_nonsexchr, file = "/scratch/gtwa/Day/sex_prediction_model/Rn7_Xist_add/processed_data/lr_sexchr_comparison/nonsexchr_only_lr.RDS")


# Predictions ####
# all features model
lr_all_prediction <- predict(lr_model_all, testing[,1:(ncol(testing)-3)], type = "response")
# sex chromosome features only
lr_sexchr_prediction <- predict(lr_model_sexchr, testing[,1:(ncol(testing)-3)], type = "response")
# non sex chromosome features only
lr_nonsexchr_prediction <- predict(lr_model_nonsexchr, testing[,1:(ncol(testing)-3)], type = "response")
## binarized classifications
lr_all_class <- ifelse(lr_all_prediction > 0.5, 1, 0) %>% as.factor()
lr_sexchr_class <- ifelse(lr_sexchr_prediction > 0.5, 1, 0) %>% as.factor()
lr_nonsexchr_class <- ifelse(lr_nonsexchr_prediction > 0.5, 1, 0) %>% as.factor()

# Create data frame for evaluation ####
eval_df <- data.frame(
  CellType = testing$CellType,
  Identity_bin = testing$Identity_bin,
  Identity_bin_char = testing$Identity_bin_char,
  lr_all_prediction = lr_all_prediction,
  lr_all_class = lr_all_class,
  lr_sexchr_prediction = lr_sexchr_prediction,
  lr_sexchr_class = lr_sexchr_class,
  lr_nonsexchr_prediction = lr_nonsexchr_prediction,
  lr_nonsexchr_class = lr_nonsexchr_class
)
## cell type specific splits
eval_df_split <- split(eval_df, f = eval_df$CellType)


# Overall Performance ####
## ROC curves ####
# calculate ROC
lr_all_roc_overall <- roc(testing$Identity_bin, lr_all_prediction, levels = c(1,0))
lr_sexchr_roc_overall <- roc(testing$Identity_bin, lr_sexchr_prediction, levels = c(1,0))
lr_nonsexchr_roc_overall <- roc(testing$Identity_bin, lr_nonsexchr_prediction, levels = c(1,0))
## save AUC values
auc_table_overall <- data.frame(AUC = c(lr_all_roc_overall[["auc"]] %>% as.numeric(),
                                        lr_sexchr_roc_overall[["auc"]] %>% as.numeric(),
                                        lr_nonsexchr_roc_overall[["auc"]] %>% as.numeric()))
rownames(auc_table_overall) <- c("LR Overall",
                                 "LR Sex Chromosome",
                                 "LR Non-Sex Chromosome")
write.csv(auc_table_overall, "/scratch/gtwa/Day/sex_prediction_model/Rn7_Xist_add/processed_data/lr_sexchr_comparison/lr_sexchr_AUC_table.csv")

# plot ROC
lr_sexchr_roccurve_overall <- ggroc(data = list(All = lr_all_roc_overall,
                                                Sexchr = lr_sexchr_roc_overall,
                                                Nonsexchr = lr_nonsexchr_roc_overall),
                                    legacy.axes = TRUE) + 
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(color = "Model") +
  theme_bw()
## save plot
ggsave(filename = "/scratch/gtwa/Day/sex_prediction_model/Rn7_Xist_add/plots/lr_sexchr_comparison/lr_sexchr_roc.pdf",
       plot = lr_sexchr_roccurve_overall,
       height = 5, width = 6)

## Confusion matrices ####
lr_all_confusion_overall <- caret::confusionMatrix(eval_df$lr_all_class, eval_df$Identity_bin)
lr_sexchr_confusion_overall <- caret::confusionMatrix(eval_df$lr_sexchr_class, eval_df$Identity_bin)
lr_nonsexchr_confusion_overall <- caret::confusionMatrix(eval_df$lr_nonsexchr_class, eval_df$Identity_bin)
# print confusion matrices to log
print("All feature Overall LR Confusion Matrix: ")
lr_all_confusion_overall
cat("\n")

print("Sex chromosome Overall LR Confusion Matrix: ")
lr_sexchr_confusion_overall
cat("\n")

print("Non-sex chromosome Overall LR Confusion Matrix: ")
lr_nonsexchr_confusion_overall
cat("\n")
# tabulate confusion matrices
confusion_list_overall <- list(
  `LR Overall` = lr_all_confusion_overall,
  `LR Sex Chromosome` = lr_sexchr_confusion_overall,
  `LR Non-Sex Chromosome` = lr_nonsexchr_confusion_overall)

confusion_df_overall <- data.frame(model = c("LR Overall", "LR Sex Chromosome", "LR Non-Sex Chromosome"),
                                   Male_correct = lapply(confusion_list_overall, function(x){x[["table"]][1,1]}) %>% unlist(),
                                   Female_incorrect = lapply(confusion_list_overall, function(x){x[["table"]][1,2]}) %>% unlist(),
                                   Male_incorrect = lapply(confusion_list_overall, function(x){x[["table"]][2,1]}) %>% unlist(),
                                   Female_correct = lapply(confusion_list_overall, function(x){x[["table"]][2,2]}) %>% unlist())

confusion_df_overall_long <- confusion_df_overall %>%
  pivot_longer(cols = Male_correct:Female_correct,
               names_to = "prediction",
               values_to = "count") %>%
  separate(col = prediction,
           into = c("reference","prediction"),
           sep = "_")

confusion_df_overall_long  <- confusion_df_overall_long  %>%
  mutate(reference = factor(reference, levels = c("Male", "Female")),
         prediction = factor(prediction, levels = c("incorrect","correct")),
         model = factor(model, levels = rev(c("LR Overall", "LR Sex Chromosome", "LR Non-Sex Chromosome"))))

# create stacked bar plot of predictions split by sex
confusion_plot_overall <- ggplot(data = confusion_df_overall_long, mapping = aes(x=model, y=count, fill=prediction)) + 
  scale_fill_manual(values = c(correct = "#00A651", incorrect = "#000000")) +
  geom_bar(position="fill", stat="identity") +
  facet_wrap( ~ reference) +
  theme_bw() +
  coord_flip() +
  labs(title = "Confusion plot", x = "Model", y = "Proportion") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
## save plot
ggsave(filename = "/scratch/gtwa/Day/sex_prediction_model/Rn7_Xist_add/plots/lr_sexchr_comparison/lr_sexchr_confusion_plot_overall.pdf",
       plot = confusion_plot_overall,
       height = 5, width = 10)

## Accuracy ####
accuracy_overall <- data.frame(Accuracy = c(lr_all_confusion_overall[["overall"]][["Accuracy"]],
                                            lr_sexchr_confusion_overall[["overall"]][["Accuracy"]],
                                            lr_nonsexchr_confusion_overall[["overall"]][["Accuracy"]]))
row.names(accuracy_overall) <- c("LR Overall", "LR Sex Chromosome", "LR Non-Sex Chromosome")
## save table
write.csv(accuracy_overall, file = "/scratch/gtwa/Day/sex_prediction_model/Rn7_Xist_add/processed_data/lr_sexchr_comparison/lr_sexchr_accuracy_overall.csv")


# CellType specific evaluation ####
## Accuracy ####
accuracy_celltype <- lapply(names(eval_df_split), function(celltype) {
  # create confusion matrices
  lr_all_confusion_celltype <- caret::confusionMatrix(eval_df_split[[celltype]]$lr_all_class, eval_df_split[[celltype]]$Identity_bin)
  lr_sexchr_confusion_celltype <- caret::confusionMatrix(eval_df_split[[celltype]]$lr_sexchr_class, eval_df_split[[celltype]]$Identity_bin)
  lr_nonsexchr_confusion_celltype <- caret::confusionMatrix(eval_df_split[[celltype]]$lr_nonsexchr_class, eval_df_split[[celltype]]$Identity_bin)
  
  # confusion matrix plots
  confusion_list <- list(`LR Overall` = lr_all_confusion_celltype,
                         `LR Sex Chromosome` = lr_sexchr_confusion_celltype,
                         `LR Non-Sex Chromosome` = lr_nonsexchr_confusion_celltype)
  confusion_df <- data.frame(model = c("LR Overall", "LR Sex Chromosome", "LR Non-Sex Chromosome"),
                             Male_correct = lapply(confusion_list, function(x){x[["table"]][1,1]}) %>% unlist(),
                             Female_incorrect = lapply(confusion_list, function(x){x[["table"]][1,2]}) %>% unlist(),
                             Male_incorrect = lapply(confusion_list, function(x){x[["table"]][2,1]}) %>% unlist(),
                             Female_correct = lapply(confusion_list, function(x){x[["table"]][2,2]}) %>% unlist())
  confusion_df_long <- confusion_df %>%
    pivot_longer(cols = Male_correct:Female_correct,
                 names_to = "prediction",
                 values_to = "count") %>%
    separate(col = prediction,
             into = c("reference","prediction"),
             sep = "_")
  confusion_df_long <- confusion_df_long  %>%
    mutate(reference = factor(reference, levels = c("Male", "Female")),
           prediction = factor(prediction, levels = c("incorrect","correct")),
           model = factor(model, levels = rev(c("LR Overall", "LR Sex Chromosome", "LR Non-Sex Chromosome"))))
  ## save table
  write.csv(confusion_df_long, file = paste0("/scratch/gtwa/Day/sex_prediction_model/Rn7_Xist_add/processed_data/lr_sexchr_comparison/lr_sexchr_confusion_table_", celltype, ".csv"))
  ## plot
  confusion_plot <- ggplot(data = confusion_df_long, mapping = aes(x=model, y=count, fill=prediction)) + 
    scale_fill_manual(values = c(correct = "#00A651", incorrect = "#000000")) +
    geom_bar(position="fill", stat="identity") +
    facet_wrap( ~ reference) +
    theme_bw() +
    coord_flip() +
    labs(title = paste0("Confuson plot: ", celltype), x = "Model", y = "Proportion") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  ## save plot
  ggsave(filename = paste0("/scratch/gtwa/Day/sex_prediction_model/Rn7_Xist_add/plots/lr_sexchr_comparison/lr_sexchr_confusion_plot_", celltype, ".pdf"),
         plot = confusion_plot,
         height = 5, width = 10)
  
  # accuracy table
  accuracy_celltype <- data.frame(Accuracy = c(lr_all_confusion_celltype[["overall"]][["Accuracy"]],
                                               lr_sexchr_confusion_celltype[["overall"]][["Accuracy"]],
                                               lr_nonsexchr_confusion_celltype[["overall"]][["Accuracy"]]))
  row.names(accuracy_celltype) <- c("LR Overall", "LR Sex Chromosome", "LR Non-Sex Chromosome")
  
  return(accuracy_celltype)
}
)
# merge tables
names(accuracy_celltype) <- names(eval_df_split)
accuracy_celltype_df <- bind_cols(accuracy_celltype)
colnames(accuracy_celltype_df) <- names(eval_df_split)
## save table
write.csv(accuracy_celltype_df, file = "/scratch/gtwa/Day/sex_prediction_model/Rn7_Xist_add/processed_data/lr_sexchr_comparison/lr_sexchr_accuracy_celltype.csv")

## ROC curves ####
auc_table_celltype <- lapply(names(eval_df_split), function(celltype) {
  # calculate ROC
  lr_all_roc_celltype <- roc(eval_df_split[[celltype]]$Identity_bin, eval_df_split[[celltype]]$lr_all_prediction, levels = c(1,0))
  lr_sexchr_roc_celltype <- roc(eval_df_split[[celltype]]$Identity_bin, eval_df_split[[celltype]]$lr_sexchr_prediction, levels = c(1,0))
  lr_nonsexchr_roc_celltype <- roc(eval_df_split[[celltype]]$Identity_bin, eval_df_split[[celltype]]$lr_nonsexchr_prediction, levels = c(1,0))
  
  # ROC curve
  roc_curve_celltype <- ggroc(data = list(All = lr_all_roc_celltype,
                                          Sexchr = lr_sexchr_roc_celltype,
                                          Nonsexchr = lr_nonsexchr_roc_celltype),
                              legacy.axes = TRUE) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    labs(color = "Model", title = celltype) +
    theme_bw()
  ## save plot
  ggsave(filename = paste0("/scratch/gtwa/Day/sex_prediction_model/Rn7_Xist_add/plots/lr_sexchr_comparison/lr_sexchr_roc_curve_", celltype, ".pdf"),
         plot = roc_curve_celltype,
         height = 5, width = 6)
  
  # AUC table
  auc_table_celltype <- data.frame(AUC = c(lr_all_roc_celltype[["auc"]] %>% as.numeric(),
                                           lr_sexchr_roc_celltype[["auc"]] %>% as.numeric(),
                                           lr_nonsexchr_roc_celltype[["auc"]] %>% as.numeric()))
  row.names(auc_table_celltype) <- c("LR Overall", "LR Sex Chromosome", "LR Non-Sex Chromosome")
  
  return(auc_table_celltype)
})
# merge tables
names(auc_table_celltype) <- names(eval_df_split)
auc_table_celltype_df <- bind_cols(auc_table_celltype)
colnames(auc_table_celltype_df) <- names(eval_df_split)
## save table
write.csv(auc_table_celltype_df, file = "/scratch/gtwa/Day/sex_prediction_model/Rn7_Xist_add/processed_data/lr_sexchr_comparison/lr_sexchr_auc_celltype.csv")

# sessionInfo ####
sessionInfo()
