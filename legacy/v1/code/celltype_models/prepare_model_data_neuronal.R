# Setup ####
## libraries
library(Seurat)
library(Boruta)
library(dplyr)

## set seed
set.seed(1234)

## load data
### Boruta object
VTA_Boruta_max2000 <- readRDS("/scratch/gtwa/Day/sex_prediction_model/data/celltype_models/VTA_Boruta_max2000_neuronal.RDS")
### training data
Rn7_VTA_training <- readRDS("/scratch/gtwa/Day/sex_prediction_model/data/celltype_models/Rn7_VTA_training_neuronal.RDS")
### training data count matrix of sex DEGs
training_count_data_sex <- readRDS("/scratch/gtwa/Day/sex_prediction_model/data/celltype_models/training_neuronal_count_subset.RDS")
### testing data
Rn7_VTA_testing <- readRDS("/scratch/gtwa/Day/sex_prediction_model/data/celltype_models/Rn7_VTA_testing_neuronal.RDS")

# Select important features ####
## Boruta information extraction
### TentativeRoughFix for tentative features
  # 54 attributes confirmed important: AC239701.1, Actb, Asb15, Atp5me, Atp5mk and 49 more;
  # 15 attributes confirmed unimportant: Abcg3l2, Calm1, Cuedc1, Eef2, ENSRNOG00000065902 and 10 more;
  # 1 tentative attributes left: Ly6h;
VTA_Boruta <- TentativeRoughFix(VTA_Boruta_max2000)
  # 55 attributes confirmed important: AC239701.1, Actb, Asb15, Atp5me, Atp5mk and 50 more;
  # 15 attributes confirmed unimportant: Abcg3l2, Calm1, Cuedc1, Eef2, ENSRNOG00000065902 and 10 more;

### extract decisions for features
VTA_Boruta_decision <- VTA_Boruta_max2000$finalDecision %>% data.frame()
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
write.csv(VTA_Boruta_decision, file = "/scratch/gtwa/Day/sex_prediction_model/data/celltype_models/Boruta_final_decisions_neuronal.csv", row.names = F)
## training and testing count data frames
write.csv(training_count_data_sex_important, file = "/scratch/gtwa/Day/sex_prediction_model/data/celltype_models/training_count_data_neuronal.csv")
write.csv(test_count_data_sex_important, file = "/scratch/gtwa/Day/sex_prediction_model/data/celltype_models/test_count_data_neuronal.csv")


# Session Info
sessionInfo()