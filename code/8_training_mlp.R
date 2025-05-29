# Set up ####
## libraries
library(RSNNS)
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
training <- read.csv("/scratch/gtwa/Day/sex_prediction_model/Rn7_Xist_add/processed_data/4_feature_selection_Boruta/training_count_data.csv", row.names = 1)
training$Identity_bin <- recode(training$Identity_bin, `0` = "Male", `1` = "Female")
training$Identity_bin <- factor(training$Identity_bin, levels = c("Male", "Female"))

# Train Model ####
# Create vectors for layer1, layer2, and decay
layer1_values <- c(1,seq(from = 5, to = 20, by = 5))
layer2_values <- c(0, round(layer1_values / 2)) %>% unique()
layer3_values <- round(layer2_values / 2)  
decay_values <- seq(from = 0, to = 0.2, by = 0.05)

# Generate all combinations
train.grid <- expand.grid(layer1 = layer1_values,
                          layer2 = layer2_values,
                          layer3 = layer3_values,
                          decay = decay_values)

# Filter rows where layer2 is smaller than layer1
train.grid <- train.grid[train.grid$layer2 < train.grid$layer1, ]

# Filter rows where layer3 is smaller than layer2
train.grid <- train.grid[train.grid$layer3 < train.grid$layer2, ]

train.ctrl <- trainControl(method = "repeatedcv",
                           number = 10,
                           repeats = 3,
                           classProbs = T,
                           allowParallel = T)
mlp_model <- caret::train(Identity_bin ~ ., data = training, 
                          method = "mlpWeightDecayML",
                          trControl = train.ctrl,
                          tuneGrid = train.grid)

mlp_model

stopCluster(cl)
# Write out files ####
saveRDS(mlp_model, file = "/scratch/gtwa/Day/sex_prediction_model/Rn7_Xist_add/processed_data/8_training_mlp/mlpwdml_model.RDS")

# print session Info ####
sessionInfo()