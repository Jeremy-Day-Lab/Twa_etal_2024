# Set up ####
## libraries
library(caret)
library(tidyverse)

## set seed
set.seed(1234)

## load data
### training dataframe
training <- read.csv("/path/to/sex_prediction_model/data/training_count_data.csv", row.names = 1)
training$Identity_bin <- as.factor(training$Identity_bin)

# Train Model ####
lr_model <- glm(Identity_bin ~ .,
                      data = training,
                      family = "binomial")
# Summary
summary(lr_model)

# Write out files ####
saveRDS(lr_model, file = "/path/to/sex_prediction_model/data/models/lr_model.RDS")

# print session Info ####
sessionInfo()
