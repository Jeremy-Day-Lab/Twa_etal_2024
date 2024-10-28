# Setup ####
## libraries
library(Seurat)
library(tidyverse)
library(caret)

## seed
set.seed(1234)

## data
Rn7_VTA <- readRDS("/path/to/sex_prediction_model/data/Rn7_VTA.RDS")

# Create split ####
## make celltype sex combination column
Rn7_VTA$CellType_Sex <- paste(Rn7_VTA$CellType,Rn7_VTA$Sex,sep = "_")

## create partition of training and test data
trainingIndex <- createDataPartition(Rn7_VTA$CellType_Sex, p = .7, 
                                     list = FALSE, 
                                     times = 1)

training.df <- Rn7_VTA@meta.data[trainingIndex,]
testing.df <- Rn7_VTA@meta.data[-trainingIndex,]

training_cells <- rownames(training.df)
testing_cells <- rownames(testing.df)

## make summary table of the number and proportion of cell type and sex 
summary.df <- table(Rn7_VTA$CellType_Sex) %>% as.data.frame()
summary.df$prop <- summary.df$Freq/sum(summary.df$Freq)

summary.df$N_cells_train <- data.frame(table(training.df$CellType_Sex))$Freq
summary.df$Proportion_train <- summary.df$N_cells_train/sum(summary.df$N_cells_train)

summary.df$N_cells_test <- data.frame(table(testing.df$CellType_Sex))$Freq
summary.df$Proportion_test <- summary.df$N_cells_test/sum(summary.df$N_cells_test)

## make subset objects
### make cell name a metadata column
Rn7_VTA@meta.data$CellName <- rownames(Rn7_VTA@meta.data)

Rn7_VTA_training <- subset(Rn7_VTA, subset = CellName %in% training_cells)
Rn7_VTA_testing <- subset(Rn7_VTA, subset = CellName %in% testing_cells)

# Save outputs ####
## save summary table
write.csv(summary.df, file = "/path/to/sex_prediction_model/data/CellType_Sex_split_summary.csv")

## save cell split vectors
saveRDS(training_cells, file = "/path/to/sex_prediction_model/data/traincell_vector.RDS")
saveRDS(testing_cells, file = "/path/to/sex_prediction_model/data/testcell_vector.RDS")

## save Rn7 training and testing objects 
saveRDS(Rn7_VTA_training, file = "/path/to/sex_prediction_model/data/Rn7_VTA_training.RDS")
saveRDS(Rn7_VTA_testing, file = "/path/to/sex_prediction_model/data/Rn7_VTA_testing.RDS")

# sessionInfo ####
sessionInfo()