# Setup ####
## libraries ####
library(Seurat)
library(Boruta)
library(dplyr)

## set seed ####
set.seed(1234)

## data ####
# training counts
training <- readRDS("/scratch/gtwa/Day/sex_prediction_model/Rn7_Xist_add/processed_data/12_create_splits_celltype_specific/training_neuronal_count_subset.RDS")
# testing counts
testing <- readRDS("/scratch/gtwa/Day/sex_prediction_model/Rn7_Xist_add/processed_data/12_create_splits_celltype_specific/Rn7_VTA_testing_neuronal.RDS")

# Feature Selection ####
VTA_Boruta <- Boruta(Identity_bin ~ ., data = training, maxRuns = 3000, doTrace = 2)
# Save output
saveRDS(VTA_Boruta,"/scratch/gtwa/Day/sex_prediction_model/Rn7_Xist_add/processed_data/neuronal/4_feature_selection_Boruta/VTA_Boruta_max3000.RDS")

# Post Selection Processing ####
# print summary
VTA_Boruta

# perform rough fix 
VTA_Boruta <- TentativeRoughFix(VTA_Boruta)
VTA_Boruta

# extract decisions for features
VTA_Boruta_decision <- VTA_Boruta$finalDecision %>% data.frame()
colnames(VTA_Boruta_decision) <- c("finalDecision")
VTA_Boruta_decision$variable <- rownames(VTA_Boruta_decision)
VTA_Boruta_decision <- VTA_Boruta_decision %>% select(variable, finalDecision) # reorder for my sanity

# Get column names with confirmed values
confirmed_columns <- VTA_Boruta_decision$variable[VTA_Boruta_decision$finalDecision == "Confirmed"]

# Create count data subsets ####
## Create training count data subsets
training_important <- training[, colnames(training) %in% confirmed_columns]
### Add identity bin back in
training_important$Identity_bin <- training$Identity_bin

## Create testing count data subset
count_data <- as.data.frame(t(as.matrix(GetAssayData(object = testing,slot = "data",assay = "RNA"))))
### check if rownames are in same order as the identities vector 
all(names(Idents(testing)) == rownames(count_data)) #TRUE 
### Now just create a new column in the count_data_full for identity
count_data$Identity <- as.factor(testing$Sex)
### Convert identitity to binary code. 
count_data$Identity_bin <- ifelse(count_data$Identity  == "Female",
                                  1, #Females are 1
                                  0) #Males are 0
### Create a subset of count data containing the same features as the training set
test_important <- count_data[,colnames(count_data) %in% colnames(training_important)]

# Save outputs ####
## Boruta variable final decisions
write.csv(VTA_Boruta_decision, file = "/scratch/gtwa/Day/sex_prediction_model/Rn7_Xist_add/processed_data/neuronal/4_feature_selection_Boruta/Boruta_final_decisions.csv", row.names = F)
## training and testing count data frames
write.csv(training_important, file = "/scratch/gtwa/Day/sex_prediction_model/Rn7_Xist_add/processed_data/neuronal/4_feature_selection_Boruta/training_count_data.csv")
write.csv(test_important, file = "/scratch/gtwa/Day/sex_prediction_model/Rn7_Xist_add/processed_data/neuronal/4_feature_selection_Boruta/test_count_data.csv")

# sessionInfo ####
sessionInfo()

# Session Info ####
sessionInfo()
