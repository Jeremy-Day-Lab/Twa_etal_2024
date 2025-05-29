# Set up ####
## libraries
library(caret)
library(tidyverse)

## set seed
set.seed(1234)

## load data
### training dataframe
training <- read.csv("/scratch/gtwa/Day/sex_prediction_model/data/celltype_models/training_count_data_neuronal.csv", row.names = 1)
training$Identity_bin <- as.factor(training$Identity_bin)

# Train Model ####
lr_model <- glm(Identity_bin ~ .,
                data = training,
                family = "binomial")
# Summary
summary(lr_model)

# Write out files ####
saveRDS(lr_model, file = "/scratch/gtwa/Day/sex_prediction_model/data/celltype_models/models/lr_model_neuronal.RDS")

# print session Info ####
sessionInfo()