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
## chrY expression based prediction
chrY_count <- y.count(testing_obj)
## Xist expression based prediction
Xist_count <- Xist.count(testing_obj)

## combined predictions
XistChrY_combined <- data.frame(
  "ChrY_count" = chrY_count,
  "Xist_count" = Xist_count,
  "ChrY_class" = ifelse(chrY_count > 0, "Male", "Female"),
  "Xist_class" = ifelse(Xist_count > 0, "Female", "Male"),
  "true_class" = testing$Identity_bin_char,
  "CellType" = testing$CellType
)
XistChrY_combined <- XistChrY_combined %>% 
  mutate(CellType = case_when(CellType %in% c("Glut-Neuron-1", "Glut-Neuron-2", "Glut-Neuron-3", "GABA-Neuron-1", "GABA-Neuron-2", "GABA-Neuron-3", "DA-Neuron") ~ "Neuronal",
                              CellType %in% c("Olig-1", "Olig-2", "Olig-3", "Astrocyte", "Polydendrocyte", "Microglia", "OPC-Olig-1", "Mural", "Endothelial") ~ "Glial"))


# Evaluate ####
## Overall ####
### How often do Xist and Chr Y based classification agree?
XistChrY_combined <- XistChrY_combined %>% 
  mutate(agreement = ifelse(ChrY_class == Xist_class, TRUE, FALSE))
agreement_prop <- (sum(XistChrY_combined$agreement) / length(XistChrY_combined$agreement)) * 100
print(paste0("Proportion of classifications in agreement: ", agreement_prop, "%"))
print(table(XistChrY_combined$agreement))


### In cases where both classifiers agree, what is the accuracy of those classifications?
XistChrY_combined <- XistChrY_combined %>% 
  mutate(correct = case_when(agreement == TRUE ~ ChrY_class == true_class,
                             agreement == FALSE ~ NA))
agreement_accuracy <- (sum(XistChrY_combined %>% filter(agreement == TRUE) %>% pull(correct)) / sum(XistChrY_combined$agreement)) * 100
print(paste0("Accuracy of classifications in agreement: ", agreement_accuracy, "%"))
print(table(XistChrY_combined %>% filter(agreement == TRUE) %>% pull(correct)))

## Neuronal ####
XistChrY_combined_neuronal <- XistChrY_combined %>% 
  filter(CellType == "Neuronal")

agreement_prop_neuronal <- (sum(XistChrY_combined_neuronal$agreement) / length(XistChrY_combined_neuronal$agreement)) * 100
print(paste0("Proportion of neuronal classifications in agreement: ", agreement_prop_neuronal, "%"))
print(table(XistChrY_combined_neuronal$agreement))

agreement_accuracy_neuronal <- (sum(XistChrY_combined_neuronal %>% filter(agreement == TRUE) %>% pull(correct)) / sum(XistChrY_combined_neuronal$agreement)) * 100
print(paste0("Accuracy of classifications in agreement: ", agreement_accuracy_neuronal, "%"))
print(table(XistChrY_combined_neuronal %>% filter(agreement == TRUE) %>% pull(correct)))

## Non-Neuronal ####
XistChrY_combined_nonneuronal <- XistChrY_combined %>% 
  filter(CellType != "Neuronal")

agreement_prop_nonneuronal <- (sum(XistChrY_combined_nonneuronal$agreement) / length(XistChrY_combined_nonneuronal$agreement)) * 100
print(paste0("Proportion of neuronal classifications in agreement: ", agreement_prop_nonneuronal, "%"))
print(table(XistChrY_combined_nonneuronal$agreement))

agreement_accuracy_nonneuronal <- (sum(XistChrY_combined_nonneuronal %>% filter(agreement == TRUE) %>% pull(correct)) / sum(XistChrY_combined_nonneuronal$agreement)) * 100
print(paste0("Accuracy of classifications in agreement: ", agreement_accuracy_nonneuronal, "%"))
print(table(XistChrY_combined_nonneuronal %>% filter(agreement == TRUE) %>% pull(correct)))

# Coalesce Results ####
res_df <- data.frame(
  Agree = c(sum(XistChrY_combined$agreement),
            sum(XistChrY_combined_neuronal$agreement),
            sum(XistChrY_combined_nonneuronal$agreement)),
  Conflict = c(length(XistChrY_combined$agreement) - sum(XistChrY_combined$agreement),
               length(XistChrY_combined_neuronal$agreement) - sum(XistChrY_combined_neuronal$agreement),
               length(XistChrY_combined_nonneuronal$agreement) - sum(XistChrY_combined_nonneuronal$agreement)),
  ConflictProportion = c(((length(XistChrY_combined$agreement) - sum(XistChrY_combined$agreement)) / length(XistChrY_combined$agreement)) * 100,
                         ((length(XistChrY_combined_neuronal$agreement) - sum(XistChrY_combined_neuronal$agreement)) / length(XistChrY_combined_neuronal$agreement)) * 100,
                         ((length(XistChrY_combined_nonneuronal$agreement) - sum(XistChrY_combined_nonneuronal$agreement)) / length(XistChrY_combined_nonneuronal$agreement)) * 100),
  AgreeProportion = c(agreement_prop,
                 agreement_prop_neuronal,
                 agreement_prop_nonneuronal),
  AgreeAccuracy = c(agreement_accuracy,
                    agreement_accuracy_neuronal,
                    agreement_accuracy_nonneuronal)
)
row.names(res_df) <- c("Overall", "Neuronal", "NonNeuronal")
res_df

write.csv(res_df, file = "/scratch/gtwa/Day/sex_prediction_model/Rn7_Xist_add/processed_data/XistChrY_combined/XistChrY_combined_stats.csv")

# sessionInfo ####
sessionInfo()