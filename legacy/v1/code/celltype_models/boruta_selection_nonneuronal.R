# Setup ####
## libraries
library(Boruta)
## seed
set.seed(1234)
## load data
count_data_subset <- readRDS("/scratch/gtwa/Day/sex_prediction_model/data/celltype_models/training_nonneuronal_count_subset.RDS")

# Feature selection ####
VTA_Boruta <- Boruta(Identity_bin ~ ., data = count_data_subset, maxRuns = 2000, doTrace = 2)

# Save output ####
saveRDS(VTA_Boruta,"/scratch/gtwa/Day/sex_prediction_model/data/celltype_models/VTA_Boruta_max2000_nonneuronal.RDS")

# Session info ####
sessionInfo()