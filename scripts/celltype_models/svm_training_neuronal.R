# Setup ####
library(kernlab)
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
### training data frame
training <- read.csv("/scratch/gtwa/Day/sex_prediction_model/data/celltype_models/training_count_data_neuronal.csv", row.names = 1)
training$Identity_bin <- recode(training$Identity_bin, `0` = "Male", `1` = "Female")
training$Identity_bin <- factor(training$Identity_bin, levels = c("Male", "Female"))

# Train Model ####
## Create vectors for sigma and C
sigma_values <- c(0.0001,0.001,0.01,0.1,1)
C_values <- c(0.01,0.1,0.2,0.5,1,1.5,2,5,10)

## Generate training grid
train.grid <- expand.grid(sigma = sigma_values, C = C_values)

## create training control parameters
train.ctrl <- trainControl(method = "repeatedcv",
                           number = 10,
                           repeats = 3,
                           classProbs = T,
                           allowParallel = T)

svm_model <- train(Identity_bin ~ ., data = training, 
                   method = "svmRadial",
                   verbose = T,
                   trControl = train.ctrl,
                   tuneGrid = train.grid)

svm_model

stopCluster(cl)
# Write out files ####
saveRDS(svm_model, file = "/scratch/gtwa/Day/sex_prediction_model/data/celltype_models/models/svm_rad_model_neuronal.RDS")

# print session Info ####
sessionInfo()