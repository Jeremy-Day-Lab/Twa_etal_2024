# Setup ####
## libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyr)
library(infotheo)
library(mpmi)

## set seed
set.seed(1234)

## data
### VTA object
Rn7_VTA <- readRDS("/path/to/sex_prediction_model/data/Rn7_VTA.RDS")
#### training object
Rn7_VTA.training <- readRDS("/path/to/sex_prediction_model/data/Rn7_VTA_training.RDS")
#### testing object
Rn7_VTA.testing <- readRDS("/path/to/sex_prediction_model/data/Rn7_VTA_testing.RDS")

## genes
selected.genes <- read.csv("/path/to/sex_prediction_model/data/Boruta_final_decisions.csv") %>% 
  filter(finalDecision == "Confirmed") %>% 
  pull("variable")

# Pre-processing data ####
## global
Rn7_VTA@meta.data <- Rn7_VTA@meta.data %>%
  mutate(Neuronal = case_when(CellType %in% c("Glut-Neuron-1", "Glut-Neuron-2", "Glut-Neuron-3", "GABA-Neuron-1", "GABA-Neuron-2", "GABA-Neuron-3", "DA-Neuron") ~ "Neuronal",
                              CellType %in% c("Olig-1", "Olig-2", "Olig-3", "Astrocyte", "Polydendrocyte", "Microglia", "OPC-Olig-1", "Mural", "Endothelial") ~ "Non-neuronal"))
Rn7_VTA@meta.data <- Rn7_VTA@meta.data %>% 
  mutate(Neuronal = factor(Neuronal, levels = c("Neuronal", "Non-neuronal")))

## training
Rn7_VTA.training@meta.data <- Rn7_VTA.training@meta.data %>%
  mutate(Neuronal = case_when(CellType %in% c("Glut-Neuron-1", "Glut-Neuron-2", "Glut-Neuron-3", "GABA-Neuron-1", "GABA-Neuron-2", "GABA-Neuron-3", "DA-Neuron") ~ "Neuronal",
                              CellType %in% c("Olig-1", "Olig-2", "Olig-3", "Astrocyte", "Polydendrocyte", "Microglia", "OPC-Olig-1", "Mural", "Endothelial") ~ "Non-neuronal"))
Rn7_VTA.training@meta.data <- Rn7_VTA.training@meta.data %>% 
  mutate(Neuronal = factor(Neuronal, levels = c("Neuronal", "Non-neuronal")))

## testing
Rn7_VTA.testing@meta.data <- Rn7_VTA.testing@meta.data %>%
  mutate(Neuronal = case_when(CellType %in% c("Glut-Neuron-1", "Glut-Neuron-2", "Glut-Neuron-3", "GABA-Neuron-1", "GABA-Neuron-2", "GABA-Neuron-3", "DA-Neuron") ~ "Neuronal",
                              CellType %in% c("Olig-1", "Olig-2", "Olig-3", "Astrocyte", "Polydendrocyte", "Microglia", "OPC-Olig-1", "Mural", "Endothelial") ~ "Non-neuronal"))
Rn7_VTA.testing@meta.data <- Rn7_VTA.testing@meta.data %>% 
  mutate(Neuronal = factor(Neuronal, levels = c("Neuronal", "Non-neuronal")))

# mpmpi calculation ####
## global ####
### extract count matrix
count.data <- as.data.frame(t(as.matrix(GetAssayData(object = Rn7_VTA,slot = "data",assay = "RNA")))) %>% 
  select(all_of(selected.genes))
### create cell lists for neuronal and non-neuronal cells
neuronal.cells <- Rn7_VTA@meta.data %>% filter(Neuronal == "Neuronal") %>% rownames()
nonneuronal.cells <- Rn7_VTA@meta.data %>% filter(Neuronal != "Neuronal") %>% rownames()

### create count dfs for neuronal and non-neuronal cells
count.data.neurons <- count.data[neuronal.cells,] 
count.data.nonneurons <- count.data[nonneuronal.cells,] 

### calculate mutual info
minfo <- lapply(selected.genes, function(gene){
  # calculate for neuronal cells
  neuronal.info <- mmi.pw(count.data.neurons %>% select(gene),
                          Rn7_VTA@meta.data[neuronal.cells,] %>% select(Sex))$mi 
  # calculate for non-neuronal cells
  nonneuronal.info <- mmi.pw(count.data.nonneurons %>% select(gene),
                             Rn7_VTA@meta.data[nonneuronal.cells,] %>% select(Sex))$mi
  # return vector of gene, minfo for neuronal, and minfo for non-neuronal
  return(c(gene, neuronal.info, nonneuronal.info))
})

### merge to dataframe
minfo.global.df <- as.data.frame(do.call(rbind,minfo))
### modify dataframe
colnames(minfo.global.df) <- c("Gene","neuronal", "nonneuronal") #set column names to be informative
minfo.global.df <- minfo.global.df %>%
  mutate(neuronal = neuronal %>% as.numeric(),
         nonneuronal = nonneuronal %>% as.numeric())
#### pivot longer
minfo.global.df.long <- minfo.global.df %>% 
  pivot_longer(cols = neuronal:nonneuronal,
               names_to = "cell.type",
               values_to = "mutual.info")


## training ####
### extract count matrix
count.data <- as.data.frame(t(as.matrix(GetAssayData(object = Rn7_VTA.training,slot = "data",assay = "RNA")))) %>% 
  select(all_of(selected.genes))
### create cell lists for neuronal and non-neuronal cells
neuronal.cells <- Rn7_VTA.training@meta.data %>% filter(Neuronal == "Neuronal") %>% rownames()
nonneuronal.cells <- Rn7_VTA.training@meta.data %>% filter(Neuronal != "Neuronal") %>% rownames()

### create count dfs for neuronal and non-neuronal cells
count.data.neurons <- count.data[neuronal.cells,] 
count.data.nonneurons <- count.data[nonneuronal.cells,] 

### calculate mutual info
minfo <- lapply(selected.genes, function(gene){
  # calculate for neuronal cells
  neuronal.info <- mmi.pw(count.data.neurons %>% select(gene),
                          Rn7_VTA.training@meta.data[neuronal.cells,] %>% select(Sex))$mi 
  # calculate for non-neuronal cells
  nonneuronal.info <- mmi.pw(count.data.nonneurons %>% select(gene),
                             Rn7_VTA.training@meta.data[nonneuronal.cells,] %>% select(Sex))$mi
  # return vector of gene, minfo for neuronal, and minfo for non-neuronal
  return(c(gene, neuronal.info, nonneuronal.info))
})

### merge to dataframe
minfo.training.df <- as.data.frame(do.call(rbind,minfo))
### modify dataframe
colnames(minfo.training.df) <- c("Gene","neuronal", "nonneuronal") #set column names to be informative
minfo.global.df <- minfo.training.df %>%
  mutate(neuronal = neuronal %>% as.numeric(),
         nonneuronal = nonneuronal %>% as.numeric())
#### pivot longer
minfo.training.df.long <- minfo.training.df %>% 
  pivot_longer(cols = neuronal:nonneuronal,
               names_to = "cell.type",
               values_to = "mutual.info")

## testing ####
### extract count matrix
count.data <- as.data.frame(t(as.matrix(GetAssayData(object = Rn7_VTA.testing,slot = "data",assay = "RNA")))) %>% 
  select(all_of(selected.genes))
### create cell lists for neuronal and non-neuronal cells
neuronal.cells <- Rn7_VTA.testing@meta.data %>% filter(Neuronal == "Neuronal") %>% rownames()
nonneuronal.cells <- Rn7_VTA.testing@meta.data %>% filter(Neuronal != "Neuronal") %>% rownames()

### create count dfs for neuronal and non-neuronal cells
count.data.neurons <- count.data[neuronal.cells,] 
count.data.nonneurons <- count.data[nonneuronal.cells,] 

### calculate mutual info
minfo <- lapply(selected.genes, function(gene){
  # calculate for neuronal cells
  neuronal.info <- mmi.pw(count.data.neurons %>% select(gene),
                          Rn7_VTA.testing@meta.data[neuronal.cells,] %>% select(Sex))$mi 
  # calculate for non-neuronal cells
  nonneuronal.info <- mmi.pw(count.data.nonneurons %>% select(gene),
                             Rn7_VTA.testing@meta.data[nonneuronal.cells,] %>% select(Sex))$mi
  # return vector of gene, minfo for neuronal, and minfo for non-neuronal
  return(c(gene, neuronal.info, nonneuronal.info))
})

### merge to dataframe
minfo.testing.df <- as.data.frame(do.call(rbind,minfo))
### modify dataframe
colnames(minfo.testing.df) <- c("Gene","neuronal", "nonneuronal") #set column names to be informative
minfo.global.df <- minfo.testing.df %>%
  mutate(neuronal = neuronal %>% as.numeric(),
         nonneuronal = nonneuronal %>% as.numeric())
#### pivot longer
minfo.testing.df.long <- minfo.testing.df %>% 
  pivot_longer(cols = neuronal:nonneuronal,
               names_to = "cell.type",
               values_to = "mutual.info")


# Write outputs ####
## global
write.csv(minfo.global.df.long, "/path/to/sex_prediction_model/data/mutualinfo.df.csv")
## training
write.csv(minfo.training.df.long, "/path/to/sex_prediction_model/data/mutualinfo.training.df.csv")
## testing
write.csv(minfo.testing.df.long, "/path/to/sex_prediction_model/data/mutualinfo.testing.df.csv")

# SessionInfo ####
sessionInfo()