# Setup ####
## libraries
library(Seurat)
library(Boruta)
library(dplyr)

## set seed
set.seed(1234)

## load data
### Boruta object
VTA_Boruta_max2000 <- readRDS("/path/to/sex_prediction_model/data/VTA_Boruta_max2000.RDS")
### training data
Rn7_VTA_training <- readRDS("/path/to/sex_prediction_model/data/Rn7_VTA_training.RDS")
### training data count matrix of sex DEGs
training_count_data_sex <- readRDS("/path/to/sex_prediction_model/data/count_data_subset_100423.RDS")
### testing data
Rn7_VTA_testing <- readRDS("/path/to/sex_prediction_model/data/Rn7_VTA_testing.RDS")

# Select important features ####
## Boruta information extraction
### TentativeRoughFix for tentative features
### Boruta ran for 1.377593 days completed its maximum number of iterations, but left 16 attributes tentative features tentative.
### we'll run `TentativeRoughFix()` to make decisions about the remaining features
VTA_Boruta <- TentativeRoughFix(VTA_Boruta_max2000)

### extract decisions for features
VTA_Boruta_decision <- VTA_Boruta$finalDecision %>% data.frame()
colnames(VTA_Boruta_decision) <- c("finalDecision")
VTA_Boruta_decision$variable <- rownames(VTA_Boruta_decision)
VTA_Boruta_decision <- VTA_Boruta_decision %>% select(variable, finalDecision) # reorder for my sanity

### Get the column names with confirmed values
confirmed_columns <- VTA_Boruta_decision$variable[VTA_Boruta_decision$finalDecision == "Confirmed"]


# Create count data subsets ####
## For training and testing data sets, we'll use dataframes of only confirmed important predictors and the outcome identity column

## Create training count data subsets
training_count_data_sex_important <- training_count_data_sex[, colnames(training_count_data_sex) %in% confirmed_columns]
### Add identity bin back in
training_count_data_sex_important$Identity_bin <- training_count_data_sex$Identity_bin

## Create testing count data subset
count_data <- as.data.frame(t(as.matrix(GetAssayData(object = Rn7_VTA_testing,slot = "data",assay = "RNA"))))
### check if rownames are in same order as the identities vector 
all(names(Idents(Rn7_VTA_testing)) == rownames(count_data)) #TRUE 
### Now just create a new column in the count_data_full for identity
count_data$Identity <- as.factor(Rn7_VTA_testing$Sex)
### Convert identitity to binary code. 
count_data$Identity_bin <- ifelse(count_data$Identity  == "Female",
                                  1, #Females are 1
                                  0) #Males are 0
### Create a subset of count data containing the same features as the training set
test_count_data_sex_important <- count_data[,colnames(count_data) %in% colnames(training_count_data_sex_important)]

# Save outputs ####
## Boruta variable final decisions
write.csv(VTA_Boruta_decision, file = "/path/to/sex_prediction_model/data/Boruta_final_decisions.csv", row.names = F)
## training and testing count data frames
write.csv(training_count_data_sex_important, file = "/path/to/sex_prediction_model/data/training_count_data.csv")
write.csv(test_count_data_sex_important, file = "/path/to/sex_prediction_model/data/test_count_data.csv")


# Session Info
sessionInfo()