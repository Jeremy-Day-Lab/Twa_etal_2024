# Setup ####
## libraries ####
library(Seurat)
library(dplyr)
library(tibble)
library(ggplot2)

## Data ####
# Original Rn7 object (does not have Xist counts)
Rn7_VTA <- readRDS("/scratch/gtwa/Day/sex_prediction_model/Rn7_Xist_add/raw_data/Rn7_VTA.RDS")

# Count matrices for Rn7+Xist alignment
F1_Xist <- Read10X("/data/project/daylab/2020-JD-0044/Rn7_Xist/F1_output_Rn7_Xist/outs/filtered_feature_bc_matrix")
F2_Xist <- Read10X("/data/project/daylab/2020-JD-0044/Rn7_Xist/F2_output_Rn7_Xist/outs/filtered_feature_bc_matrix")
M1_Xist <- Read10X("/data/project/daylab/2020-JD-0044/Rn7_Xist/M1_output_Rn7_Xist/outs/filtered_feature_bc_matrix")
M2_Xist <- Read10X("/data/project/daylab/2020-JD-0044/Rn7_Xist/M2_output_Rn7_Xist/outs/filtered_feature_bc_matrix")

# Create Xist object ####
# Create individual objects for each sample
Fem1  <- CreateSeuratObject(counts = F1_Xist,min.cells = 0,min.features = 0) #5936 nuclei, original: 5940  
Fem2  <- CreateSeuratObject(counts = F2_Xist,min.cells = 0,min.features = 0) #5274 nuclei, original: 5278 
Male1 <- CreateSeuratObject(counts = M1_Xist,min.cells = 0,min.features = 0) #3851 nuclei, original: 3841 
Male2 <- CreateSeuratObject(counts = M2_Xist,min.cells = 0,min.features = 0) #7112 nuclei, original: 7120  
## 22173 total nuclei, original: 22179

# Rename cells to match original Rn7 object
Fem1 <- RenameCells(Fem1, new.names = paste0(colnames(Fem1), "_1"))
Fem2 <- RenameCells(Fem2, new.names = paste0(colnames(Fem2), "_2"))
Male1 <- RenameCells(Male1, new.names = paste0(colnames(Male1), "_3"))
Male2 <- RenameCells(Male2, new.names = paste0(colnames(Male2), "_4"))

# Verify that new cell names match a subset of the original Rn7 object
Fem1_cells <- Cells(Fem1)[Cells(Fem1) %in% c(Rn7_VTA@meta.data %>% filter(GEM_Well == "Fem1") %>% rownames())] #5930 nuclei
Fem2_cells <- Cells(Fem2)[Cells(Fem2) %in% c(Rn7_VTA@meta.data %>% filter(GEM_Well == "Fem2") %>% rownames())] #5274 nuclei
Male1_cells <- Cells(Male1)[Cells(Male1) %in% c(Rn7_VTA@meta.data %>% filter(GEM_Well == "Male1") %>% rownames())] #3839 nuclei
Male2_cells <- Cells(Male2)[Cells(Male2) %in% c(Rn7_VTA@meta.data %>% filter(GEM_Well == "Male2") %>% rownames())] #7106 nuclei
Xist_cells <- c(Fem1_cells, Fem2_cells, Male1_cells, Male2_cells)
## 22149 total nuclei, original final object: 22170

# Extract raw Xist counts from Xist objects
Fem1_Xist <- LayerData(Fem1, "counts", cells = Fem1_cells, features = "Xist")
Fem2_Xist <- LayerData(Fem2, "counts", cells = Fem2_cells, features = "Xist")
Male1_Xist <- LayerData(Male1, "counts", cells = Male1_cells, features = "Xist")
Male2_Xist <- LayerData(Male2, "counts", cells = Male2_cells, features = "Xist")
## Merge Xist counts into a single matrix
Xist_counts <- cbind(Fem1_Xist, Fem2_Xist, Male1_Xist, Male2_Xist)

# Extract raw counts from Rn7 object
Rn7_VTA_Xist <- LayerData(Rn7_VTA, "counts", cells = Xist_cells)
## Merge Xist counts with Rn7 counts
Rn7_VTA_Xist <- rbind(Rn7_VTA_Xist, Xist_counts)

# Create Object
Rn7_VTA_Xist_obj <- CreateSeuratObject(counts = Rn7_VTA_Xist, min.cells = 0, min.features = 0)
## Re-add metadata
Rn7_VTA_Xist_obj$GEM_Well <- Rn7_VTA@meta.data[rownames(Rn7_VTA_Xist_obj@meta.data), "GEM_Well"]
Rn7_VTA_Xist_obj$Sex <- Rn7_VTA@meta.data[rownames(Rn7_VTA_Xist_obj@meta.data), "Sex"]
Rn7_VTA_Xist_obj$CellType <-  Rn7_VTA@meta.data[rownames(Rn7_VTA_Xist_obj@meta.data), "CellType"]
## Normalize data
Rn7_VTA_Xist_obj <- NormalizeData(Rn7_VTA_Xist_obj, normalization.method = "LogNormalize", scale.factor = 10000)
## Scale data
Rn7_VTA_Xist_obj <- ScaleData(Rn7_VTA_Xist_obj, verbose = FALSE)

# Save object ####
saveRDS(Rn7_VTA_Xist_obj, file = "/scratch/gtwa/Day/sex_prediction_model/Rn7_Xist_add/processed_data/Rn7_Xist_VTA.RDS")

# sessionInfo ####
sessionInfo()