# Set up ####
## libraries
library(caret)
library(tidyverse)

## set seed
set.seed(1234)

## load data
### training dataframe
training <- read.csv("/scratch/gtwa/Day/sex_prediction_model/Rn7_Xist_add/processed_data/neuronal/4_feature_selection_Boruta/training_count_data.csv", row.names = 1)
training$Identity_bin <- as.factor(training$Identity_bin)

# Train Model ####
lr_model <- glm(Identity_bin ~ .,
                data = training,
                family = "binomial")
# Summary
summary(lr_model)

# Write out files ####
saveRDS(lr_model, file = "/scratch/gtwa/Day/sex_prediction_model/Rn7_Xist_add/processed_data/neuronal/5_training_logistic_regression/lr_model.RDS")

# print session Info ####
sessionInfo()