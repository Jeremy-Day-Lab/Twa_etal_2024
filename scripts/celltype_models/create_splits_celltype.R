# Setup ####
## libraries
library(Seurat)
library(tidyverse)

## seed
set.seed(1234)

## data
### training
Rn7_VTA.training <- readRDS("/scratch/gtwa/Day/sex_prediction_model/data/Rn7_VTA_training.RDS")
### testing
Rn7_VTA.testing <- readRDS("/scratch/gtwa/Day/sex_prediction_model/data/Rn7_VTA_testing.RDS")
### previous DEG results
Sex_DEGs <- read.csv("/scratch/gtwa/Day/sex_prediction_model/data/DEGs_sex_cluster_VTA_sig.csv")
### GTF
# Rn7_gtf <- rtracklayer::import("/scratch/gtwa/Day/sex_prediction_model/data/Rattus_norvegicus.mRatBN7.2.105.gtf")
# Rn7_gtf <- as.data.frame(Rn7_gtf)


# pre-processing ####
## make neuronal and non-neuronal cell labels
Rn7_VTA.training@meta.data <- Rn7_VTA.training@meta.data %>%
  mutate(Neuronal = case_when(CellType %in% c("Glut-Neuron-1", "Glut-Neuron-2", "Glut-Neuron-3", "GABA-Neuron-1", "GABA-Neuron-2", "GABA-Neuron-3", "DA-Neuron") ~ "Neuronal",
                              CellType %in% c("Olig-1", "Olig-2", "Olig-3", "Astrocyte", "Polydendrocyte", "Microglia", "OPC-Olig-1", "Mural", "Endothelial") ~ "Non-neuronal")) # make a new column to designate neuronal and non-neuronal
Rn7_VTA.training@meta.data <- Rn7_VTA.training@meta.data %>% 
  mutate(Neuronal = factor(Neuronal, levels = c("Neuronal", "Non-neuronal"))) # make cell type and neuronal/non-neuronal factors 

Rn7_VTA.testing@meta.data <- Rn7_VTA.testing@meta.data %>%
  mutate(Neuronal = case_when(CellType %in% c("Glut-Neuron-1", "Glut-Neuron-2", "Glut-Neuron-3", "GABA-Neuron-1", "GABA-Neuron-2", "GABA-Neuron-3", "DA-Neuron") ~ "Neuronal",
                              CellType %in% c("Olig-1", "Olig-2", "Olig-3", "Astrocyte", "Polydendrocyte", "Microglia", "OPC-Olig-1", "Mural", "Endothelial") ~ "Non-neuronal")) # make a new column to designate neuronal and non-neuronal
Rn7_VTA.testing@meta.data <- Rn7_VTA.testing@meta.data %>% 
  mutate(Neuronal = factor(Neuronal, levels = c("Neuronal", "Non-neuronal"))) # make cell type and neuronal/non-neuronal factors 


# split DEG list by neuronal and non-neuronal cluster results ####
### neuronal
neuronal.degs <- Sex_DEGs %>%
  filter(Cluster %in% c("Glut-Neuron-1", "Glut-Neuron-2", "Glut-Neuron-3", "GABA-Neuron-1", "GABA-Neuron-2", "GABA-Neuron-3", "DA-Neuron")) %>% 
  pull(GeneName) %>% 
  unique() # 70 genes

### non-neuronal
nonneuronal.degs <- Sex_DEGs %>% filter(Cluster %in% c("Olig-1", "Olig-2", "Olig-3", "Astrocyte", "Polydendrocyte", "Microglia", "OPC-Olig-1", "Mural", "Endothelial")) %>% 
  pull(GeneName) %>% 
  unique() # 364 genes


# make subsets for cell type ####
## training
Rn7_VTA.training@active.ident <- Rn7_VTA.training$Neuronal
Rn7_VTA.training.neuronal <- subset(x = Rn7_VTA.training, idents = "Neuronal")
Rn7_VTA.training.nonneuronal <- subset(x = Rn7_VTA.training, idents = "Non-neuronal")

## testing
Rn7_VTA.testing@active.ident <- Rn7_VTA.testing$Neuronal
Rn7_VTA.testing.neuronal <- subset(x = Rn7_VTA.testing, idents = "Neuronal")
Rn7_VTA.testing.nonneuronal <- subset(x = Rn7_VTA.testing, idents = "Non-neuronal")

# make subsets of count data ####
## training 
### neuronal
#### First pull the count data
training.neuronal.count_data_Full <- as.data.frame(t(as.matrix(GetAssayData(object = Rn7_VTA.training.neuronal,slot = "data",assay = "RNA"))))
#### check if rownames are in same order as the identities vector 
all(names(Idents(Rn7_VTA.training.neuronal)) == rownames(training.neuronal.count_data_Full)) #TRUE 
#### Now just create a new column in the training.neuronal.count_data_Full for identity
training.neuronal.count_data_Full$Identity <- as.factor(Rn7_VTA.training.neuronal$Sex)
#### Convert identitity to binary code. 
training.neuronal.count_data_Full$Identity_bin <- ifelse(training.neuronal.count_data_Full$Identity  == "Female",
                                                         1, #Females are 1
                                                         0) #Males are 0
#### Create a subset of count data full containing all of the Sex_DEGs
training.neuronal.count_data_subset <- training.neuronal.count_data_Full[,c(neuronal.degs,"Identity_bin")]


### non-neuronal
#### First pull the count data
training.nonneuronal.count_data_Full <- as.data.frame(t(as.matrix(GetAssayData(object = Rn7_VTA.training.nonneuronal,slot = "data",assay = "RNA"))))
#### check if rownames are in same order as the identities vector 
all(names(Idents(Rn7_VTA.training.nonneuronal)) == rownames(training.nonneuronal.count_data_Full)) #TRUE 
#### Now just create a new column in the training.nonneuronal.count_data_Full for identity
training.nonneuronal.count_data_Full$Identity <- as.factor(Rn7_VTA.training.nonneuronal$Sex)
#### Convert identitity to binary code. 
training.nonneuronal.count_data_Full$Identity_bin <- ifelse(training.nonneuronal.count_data_Full$Identity  == "Female",
                                                            1, #Females are 1
                                                            0) #Males are 0
#### Create a subset of count data full containing all of the Sex_DEGs
training.nonneuronal.count_data_subset <- training.nonneuronal.count_data_Full[,c(nonneuronal.degs,"Identity_bin")]


# save outputs ####
## seurat objects
saveRDS(Rn7_VTA.training.neuronal, "/scratch/gtwa/Day/sex_prediction_model/data/celltype_models/Rn7_VTA_training_neuronal.RDS")
saveRDS(Rn7_VTA.training.nonneuronal, "/scratch/gtwa/Day/sex_prediction_model/data/celltype_models/Rn7_VTA_training_nonneuronal.RDS")

saveRDS(Rn7_VTA.testing.neuronal, "/scratch/gtwa/Day/sex_prediction_model/data/celltype_models/Rn7_VTA_testing_neuronal.RDS")
saveRDS(Rn7_VTA.testing.nonneuronal, "/scratch/gtwa/Day/sex_prediction_model/data/celltype_models/Rn7_VTA_testing_nonneuronal.RDS")

## count data subsets
saveRDS(training.neuronal.count_data_subset, "/scratch/gtwa/Day/sex_prediction_model/data/celltype_models/training_neuronal_count_subset.RDS")
saveRDS(training.nonneuronal.count_data_subset, "/scratch/gtwa/Day/sex_prediction_model/data/celltype_models/training_nonneuronal_count_subset.RDS")

# sessionInfo ####
sessionInfo()