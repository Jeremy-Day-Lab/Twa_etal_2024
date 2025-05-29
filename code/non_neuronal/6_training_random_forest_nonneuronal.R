# Setup ####
library(randomForest)
library(ROCR)
library(caret)
library(tidyverse)
library(doParallel)

## set seed
set.seed(1234)
## setup parallel
cl <- makePSOCKcluster(10)
registerDoParallel(cl)

## load data
### training dataframe
training <- read.csv("/scratch/gtwa/Day/sex_prediction_model/Rn7_Xist_add/processed_data/non_neuronal/4_feature_selection_Boruta/training_count_data.csv", row.names = 1)
training$Identity_bin <- recode(training$Identity_bin, `0` = "Male", `1` = "Female")
training$Identity_bin <- factor(training$Identity_bin, levels = c("Male", "Female"))

# Train Model ####
mtry_values <- seq(from = 2, to = 334, by = 12)

# Generate all combinations
tune.grid <- expand.grid(mtry = mtry_values)

train.ctrl <- trainControl(method = "repeatedcv",
                           number = 10,
                           repeats = 3,
                           classProbs = T,
                           allowParallel = T)
rf_model <- train(Identity_bin ~ ., data = training, 
                  method = "rf",
                  trControl = train.ctrl,
                  tuneGrid = tune.grid,
                  ntree = 1000)
rf_model

stopCluster(cl)

# Write out files ####
saveRDS(rf_model, file = "/scratch/gtwa/Day/sex_prediction_model/Rn7_Xist_add/processed_data/non_neuronal/6_training_random_forest/rf_model.RDS")

# print session Info ####
sessionInfo()