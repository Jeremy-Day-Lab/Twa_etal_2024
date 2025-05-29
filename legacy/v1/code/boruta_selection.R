# Setup ####
## libraries
library(Boruta)
## seed
set.seed(1234)
## load data
count_data_subset <- readRDS("/path/to/sex_prediction_model/data/count_data_subset_100423.RDS")

# Feature selection ####
VTA_Boruta <- Boruta(Identity_bin ~ ., data = count_data_subset, maxRuns = 2000, doTrace = 2)

# Save output ####
saveRDS(VTA_Boruta,"/path/to/sex_prediction_model/data/VTA_Boruta_max2000.RDS")

# Session info ####
sessionInfo()