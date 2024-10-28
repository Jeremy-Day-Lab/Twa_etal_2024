# Setup ####
## libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(tibble)
library(tidyr)

## set seed
set.seed(1234)

## functions
### count Y transcripts
y.count <- function(srat.object){
  # load gene expression count data
  srat.counts <- as.data.frame(t(as.matrix(GetAssayData(object = srat.object, slot = "data", assay = "RNA"))))
  
  # load gene annotations
  Rn7.gtf <- rtracklayer::import("/path/to/sex_prediction_model/data/Rattus_norvegicus.mRatBN7.2.105.gtf")
  Rn7.gtf <- as.data.frame(Rn7.gtf)
  
  # select Y chromosome genes
  y.genes <- Rn7.gtf %>% #from Rn7 gene annotations
    filter(seqnames == "Y") %>% #filter for chrY only
    select(gene_id, gene_name) %>% #select gene ids and names
    mutate(gene_name = coalesce(gene_name, gene_id)) %>% #replace NA name values with gene ID
    pull(gene_name) %>% #pull gene name values
    unique() #only unique values, returns 27 genes
  
  # create cumulative expresion of Y genes
  y.counts <- srat.counts %>%
    select(any_of(y.genes)) %>%
    rowSums() %>% data.frame()
  ### set column name
  colnames(y.counts) <- c("y.counts")
  
  return(y.counts)
}
### count Xist transcripts
Xist.count <- function(srat.object){
  # load gene expression count data
  Xist.counts <- as.data.frame(t(as.matrix(GetAssayData(object = srat.object, slot = "data", assay = "RNA")))) %>%
    pull("ENSRNOG00000065796") # pull out Xist (named as ENSRNOG00000065796)
  
  return(Xist.counts)
}

## data
### VTA object
Rn7_VTA <- readRDS("/path/to/sex_prediction_model/data/Rn7_VTA.RDS")
#### training object
Rn7_VTA.training <- readRDS("/path/to/sex_prediction_model/data/Rn7_VTA_training.RDS")
#### testing object
Rn7_VTA.testing <- readRDS("/path/to/sex_prediction_model/data/Rn7_VTA_testing.RDS")

## genes
### Rn7 gene annotations
Rn7.gtf <- rtracklayer::import("/path/to/sex_prediction_model/data/Rattus_norvegicus.mRatBN7.2.105.gtf")
Rn7.gtf <- as.data.frame(Rn7.gtf)
### selected genes
selected.genes <- read.csv("/path/to/sex_prediction_model/data/Boruta_final_decisions.csv") %>% 
  filter(finalDecision == "Confirmed") %>% 
  pull("variable")
### Y chromosome genes
y.genes <- Rn7.gtf %>% #from Rn7 gene annotations
  filter(seqnames == "Y") %>% #filter for chrY only
  select(gene_id, gene_name) %>% #select gene ids and names
  mutate(gene_name = coalesce(gene_name, gene_id)) %>% #replace NA name values with gene ID
  pull(gene_name) %>% #pull gene name values
  unique() #only unique values, returns 27 genes


# Pre-processing data ####
## count Y genes
Rn7_VTA$y.counts <- y.count(Rn7_VTA) # whole object
Rn7_VTA.training$y.counts <- y.count(Rn7_VTA.training) # training object
Rn7_VTA.testing$y.counts <- y.count(Rn7_VTA.testing) # testing object
## count X genes
Rn7_VTA$Xist.counts <- Xist.count(Rn7_VTA)
Rn7_VTA.training$Xist.counts <- Xist.count(Rn7_VTA.training)
Rn7_VTA.testing$Xist.counts <- Xist.count(Rn7_VTA.testing)

## Make neuronal/non-neuronal cell type column
Rn7_VTA@meta.data <- Rn7_VTA@meta.data %>%
  mutate(Neuronal = case_when(CellType %in% c("Glut-Neuron-1", "Glut-Neuron-2", "Glut-Neuron-3", "GABA-Neuron-1", "GABA-Neuron-2", "GABA-Neuron-3", "DA-Neuron") ~ "Neuronal",
                              CellType %in% c("Olig-1", "Olig-2", "Olig-3", "Astrocyte", "Polydendrocyte", "Microglia", "OPC-Olig-1", "Mural", "Endothelial") ~ "Non-neuronal")) # make a new column to designate neuronal and non-neuronal
Rn7_VTA@meta.data <- Rn7_VTA@meta.data %>% 
  mutate(Neuronal = factor(Neuronal, levels = c("Neuronal", "Non-neuronal")),
         CellType = factor(CellType, levels = c("DA-Neuron","Glut-Neuron-1", "Glut-Neuron-2", "Glut-Neuron-3", "GABA-Neuron-1", "GABA-Neuron-2", "GABA-Neuron-3",
                                                "Olig-1", "Olig-2", "Olig-3", "Astrocyte","Microglia","Polydendrocyte","OPC-Olig-1", "Mural", "Endothelial"))) # make cell type and neuronal/non-neuronal factors 

Rn7_VTA.training@meta.data <- Rn7_VTA.training@meta.data %>%
  mutate(Neuronal = case_when(CellType %in% c("Glut-Neuron-1", "Glut-Neuron-2", "Glut-Neuron-3", "GABA-Neuron-1", "GABA-Neuron-2", "GABA-Neuron-3", "DA-Neuron") ~ "Neuronal",
                              CellType %in% c("Olig-1", "Olig-2", "Olig-3", "Astrocyte", "Polydendrocyte", "Microglia", "OPC-Olig-1", "Mural", "Endothelial") ~ "Non-neuronal")) # make a new column to designate neuronal and non-neuronal
Rn7_VTA.training@meta.data <- Rn7_VTA.training@meta.data %>% 
  mutate(Neuronal = factor(Neuronal, levels = c("Neuronal", "Non-neuronal")),
         CellType = factor(CellType, levels = c("DA-Neuron","Glut-Neuron-1", "Glut-Neuron-2", "Glut-Neuron-3", "GABA-Neuron-1", "GABA-Neuron-2", "GABA-Neuron-3",
                                                "Olig-1", "Olig-2", "Olig-3", "Astrocyte","Microglia","Polydendrocyte","OPC-Olig-1", "Mural", "Endothelial"))) # make cell type and neuronal/non-neuronal factors 

Rn7_VTA.testing@meta.data <- Rn7_VTA.testing@meta.data %>%
  mutate(Neuronal = case_when(CellType %in% c("Glut-Neuron-1", "Glut-Neuron-2", "Glut-Neuron-3", "GABA-Neuron-1", "GABA-Neuron-2", "GABA-Neuron-3", "DA-Neuron") ~ "Neuronal",
                              CellType %in% c("Olig-1", "Olig-2", "Olig-3", "Astrocyte", "Polydendrocyte", "Microglia", "OPC-Olig-1", "Mural", "Endothelial") ~ "Non-neuronal")) # make a new column to designate neuronal and non-neuronal
Rn7_VTA.testing@meta.data <- Rn7_VTA.testing@meta.data %>% 
  mutate(Neuronal = factor(Neuronal, levels = c("Neuronal", "Non-neuronal")),
         CellType = factor(CellType, levels = c("DA-Neuron","Glut-Neuron-1", "Glut-Neuron-2", "Glut-Neuron-3", "GABA-Neuron-1", "GABA-Neuron-2", "GABA-Neuron-3",
                                                "Olig-1", "Olig-2", "Olig-3", "Astrocyte","Microglia","Polydendrocyte","OPC-Olig-1", "Mural", "Endothelial"))) # make cell type and neuronal/non-neuronal factors 

## Sex distributions
### all cells
all.sex.distb.df.long <- table(Rn7_VTA$CellType, Rn7_VTA$Sex) %>% # make a table of cell counts per cell type by sex
  as.data.frame.matrix() %>% # convert to dataframe
  rownames_to_column("CellType") %>% # make cell type a column
  mutate(CellType = factor(CellType, levels = c("DA-Neuron","Glut-Neuron-1", "Glut-Neuron-2", "Glut-Neuron-3", "GABA-Neuron-1", "GABA-Neuron-2", "GABA-Neuron-3",
                                                "Olig-1", "Olig-2", "Olig-3", "Astrocyte","Microglia","Polydendrocyte","OPC-Olig-1", "Mural", "Endothelial"))) %>% # make cell type a factor
  pivot_longer(cols = Female:Male, names_to = "Sex", values_to = "count") # pivot longer for plotting
### training cells
training.sex.distb.df.long <- table(Rn7_VTA.training$CellType, Rn7_VTA.training$Sex) %>% # make a table of cell counts per cell type by sex
  as.data.frame.matrix() %>% # convert to dataframe
  rownames_to_column("CellType") %>% # make cell type a column
  mutate(CellType = factor(CellType, levels = c("DA-Neuron","Glut-Neuron-1", "Glut-Neuron-2", "Glut-Neuron-3", "GABA-Neuron-1", "GABA-Neuron-2", "GABA-Neuron-3",
                                                "Olig-1", "Olig-2", "Olig-3", "Astrocyte","Microglia","Polydendrocyte","OPC-Olig-1", "Mural", "Endothelial"))) %>% # make cell type a factor
  pivot_longer(cols = Female:Male, names_to = "Sex", values_to = "count") # pivot longer for plotting
### testing cells
testing.sex.distb.df.long <- table(Rn7_VTA.testing$CellType, Rn7_VTA.testing$Sex) %>% # make a table of cell counts per cell type by sex
  as.data.frame.matrix() %>% # convert to dataframe
  rownames_to_column("CellType") %>% # make cell type a column
  mutate(CellType = factor(CellType, levels = c("DA-Neuron","Glut-Neuron-1", "Glut-Neuron-2", "Glut-Neuron-3", "GABA-Neuron-1", "GABA-Neuron-2", "GABA-Neuron-3",
                                                "Olig-1", "Olig-2", "Olig-3", "Astrocyte","Microglia","Polydendrocyte","OPC-Olig-1", "Mural", "Endothelial"))) %>% # make cell type a factor
  pivot_longer(cols = Female:Male, names_to = "Sex", values_to = "count") # pivot longer for plotting

## ncells per cluster
### all cells
all.ncells.df <- table(Rn7_VTA$CellType) %>% # make a table of cell type counts
  data.frame() %>% # convert to dataframe
  rename("CellType" = "Var1",
         "Cells" = "Freq") %>% # rename columns
  mutate(CellType = factor(CellType, levels = c("DA-Neuron","Glut-Neuron-1", "Glut-Neuron-2", "Glut-Neuron-3", "GABA-Neuron-1", "GABA-Neuron-2", "GABA-Neuron-3",
                                                "Olig-1", "Olig-2", "Olig-3", "Astrocyte","Microglia","Polydendrocyte","OPC-Olig-1", "Mural", "Endothelial"))) # make cell type a factor
### training cells
training.ncells.df <- table(Rn7_VTA.training$CellType) %>% # make a table of cell type counts
  data.frame() %>% # convert to dataframe
  rename("CellType" = "Var1",
         "Cells" = "Freq") %>% # rename columns
  mutate(CellType = factor(CellType, levels = c("DA-Neuron","Glut-Neuron-1", "Glut-Neuron-2", "Glut-Neuron-3", "GABA-Neuron-1", "GABA-Neuron-2", "GABA-Neuron-3",
                                                "Olig-1", "Olig-2", "Olig-3", "Astrocyte","Microglia","Polydendrocyte","OPC-Olig-1", "Mural", "Endothelial"))) # make cell type a factor
### testing cells
testing.ncells.df <- table(Rn7_VTA.testing$CellType) %>% # make a table of cell type counts
  data.frame() %>% # convert to dataframe
  rename("CellType" = "Var1",
         "Cells" = "Freq") %>% # rename columns
  mutate(CellType = factor(CellType, levels = c("DA-Neuron","Glut-Neuron-1", "Glut-Neuron-2", "Glut-Neuron-3", "GABA-Neuron-1", "GABA-Neuron-2", "GABA-Neuron-3",
                                                "Olig-1", "Olig-2", "Olig-3", "Astrocyte","Microglia","Polydendrocyte","OPC-Olig-1", "Mural", "Endothelial"))) # make cell type a factor


# Plots ####
## whole object
### clustering
#### UMAP
Rn7_VTA@active.ident <- Rn7_VTA$CellType
all.dimplot.clusters <- DimPlot(Rn7_VTA)
#### Number of cells per cluster
all.ncells.clusters <- ggplot(data = all.ncells.df, mapping = aes(x = CellType, y = Cells, fill = CellType)) +
  geom_bar(stat = "identity") +
  scale_y_log10() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 15),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 15),
        panel.grid = element_blank()) +
  NoLegend()

### distribution of sexes
#### UMAP
Rn7_VTA@active.ident <- factor(Rn7_VTA$Sex, levels = c("Female", "Male")) # set active identity to sex for this plot
all.dimplot.sexes <- DimPlot(Rn7_VTA, shuffle = T)
#### proportions per cluster
all.prop.sexes <- ggplot(data = all.sex.distb.df.long, mapping = aes(x = CellType, y = count, fill = Sex)) +
  geom_bar(position="fill", stat="identity") +
  labs(x = "Cell Type", y = "Proportion") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 15),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 15),
        panel.grid = element_blank()) +
  NoLegend()

Rn7_VTA@active.ident <- Rn7_VTA$CellType # set active identity back to cell type

### Y gene counts
#### feature plot
all.dim.plot.y <- FeaturePlot(Rn7_VTA,
                              features = "y.counts",
                              split.by = "Sex",
                              cols = c("lightgrey","blue"),
                              keep.scale = "all")
#### violin
all.vln.plot.y <- VlnPlot(Rn7_VTA,
                          features = "y.counts",
                          split.by = "Sex",
                          pt.size = 0) +
  labs(x = "Cell Type", y = "ChrY Gene Expression") +
  NoLegend()

### Xist
#### feature plot
all.dim.plot.xist <- FeaturePlot(Rn7_VTA,
                                 features = "Xist.counts",
                                 split.by = "Sex",
                                 cols = c("lightgrey","red"),
                                 keep.scale = "all") 

#### violin
all.vln.plot.xist <- VlnPlot(Rn7_VTA,
                             features = "Xist.counts",
                             split.by = "Sex",
                             pt.size = 0) +
  labs(x = "Cell Type", y = "Xist Gene Expression") +
  NoLegend()

### DEGs
all.vln.plot.degs <- VlnPlot(Rn7_VTA,
                             features = c("ENSRNOG00000060617",
                                          "ENSRNOG00000065796",
                                          "Kdm5d",
                                          "Eif2s3y",
                                          "Kdm6a",
                                          "Cntnap2",
                                          "Crip1"),
                             split.by = "Sex",
                             pt.size = 0,
                             stack = T,
                             flip = T) +
  labs(x = "Cell Type") +
  theme(legend.position = "top")

### gene counts
#### distribution
all.gene.count.distb <- ggplot(data = Rn7_VTA@meta.data, mapping = aes(x = nFeature_RNA, fill = Neuronal)) +
  geom_density(alpha = 0.7) +
  labs(x = "Number of Genes", y = "Density") +
  theme_bw() +
  theme(text = element_text(size = 15),
        axis.text = element_text(size = 10),
        legend.position = c(0.9, 0.9),
        legend.title= element_blank()) +
  NoGrid()


## training object
### clustering
#### UMAP
Rn7_VTA.training@active.ident <- Rn7_VTA.training$CellType
training.dimplot.clusters <- DimPlot(Rn7_VTA.training)
#### Number of cells per cluster
training.ncells.clusters <- ggplot(data = training.ncells.df, mapping = aes(x = CellType, y = Cells, fill = CellType)) +
  geom_bar(stat = "identity") +
  scale_y_log10() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 15),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 15),
        panel.grid = element_blank()) +
  NoLegend()

### distribution of sexes
#### UMAP
Rn7_VTA.training@active.ident <- factor(Rn7_VTA.training$Sex, levels = c("Female", "Male")) # set active ident to sex for this plot
training.dimplot.sexes <- DimPlot(Rn7_VTA.training)
#### proportions per cluster
training.prop.sexes <- ggplot(data = training.sex.distb.df.long, mapping = aes(x = CellType, y = count, fill = Sex)) +
  geom_bar(position="fill", stat="identity") +
  labs(x = "Cell Type", y = "Proportion") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 15),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 15),
        panel.grid = element_blank()) +
  NoLegend()

Rn7_VTA.training@active.ident <- Rn7_VTA.training$CellType # reset active ident to cell type

### Y gene counts
#### feature plot
training.dim.plot.y <- FeaturePlot(Rn7_VTA.training,
                                   features = "y.counts",
                                   split.by = "Sex",
                                   cols = c("lightgrey","blue"),
                                   keep.scale = "all")
#### violin
training.vln.plot.y <- VlnPlot(Rn7_VTA.training,
                               features = "y.counts",
                               split.by = "Sex",
                               pt.size = 0) +
  labs(x = "Cell Type", y = "ChrY Gene Expression") +
  NoLegend()

### Xist
#### feature plot
training.dim.plot.xist <- FeaturePlot(Rn7_VTA.training,
                                      features = "Xist.counts",
                                      split.by = "Sex",
                                      cols = c("lightgrey","red"),
                                      keep.scale = "all") 

#### violin
training.vln.plot.xist <- VlnPlot(Rn7_VTA.training,
                                  features = "Xist.counts",
                                  split.by = "Sex",
                                  pt.size = 0) +
  labs(x = "Cell Type", y = "Xist Gene Expression") +
  NoLegend()

### DEGs
training.vln.plot.degs <- VlnPlot(Rn7_VTA.training,
                                  features = c("ENSRNOG00000060617",
                                               "ENSRNOG00000065796",
                                               "Kdm5d",
                                               "Eif2s3y",
                                               "Kdm6a",
                                               "Cntnap2",
                                               "Crip1"),
                                  split.by = "Sex",
                                  pt.size = 0,
                                  stack = T,
                                  flip = T) +
  labs(x = "Cell Type") +
  theme(legend.position = "top")


### gene counts
#### distribution
training.gene.count.distb <- ggplot(data = Rn7_VTA.training@meta.data, mapping = aes(x = nFeature_RNA, fill = Neuronal)) +
  geom_density(alpha = 0.7) +
  labs(x = "Number of Genes", y = "Density") +
  theme_bw() +
  theme(text = element_text(size = 15),
        axis.text = element_text(size = 10),
        legend.position = c(0.9, 0.9),
        legend.title= element_blank()) +
  NoGrid()

## testing object
### clustering
#### UMAP
Rn7_VTA.testing@active.ident <- Rn7_VTA.testing$CellType
testing.dimplot.clusters <- DimPlot(Rn7_VTA.testing)
#### Number of cells per cluster
testing.ncells.clusters <- ggplot(data = testing.ncells.df, mapping = aes(x = CellType, y = Cells, fill = CellType)) +
  geom_bar(stat = "identity") +
  scale_y_log10() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 15),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 15),
        panel.grid = element_blank()) +
  NoLegend()

### distribution of sexes
#### UMAP
Rn7_VTA.testing@active.ident <- factor(Rn7_VTA.testing$Sex, levels = c("Female", "Male")) # set active ident to sex for this plot
testing.dimplot.sexes <- DimPlot(Rn7_VTA.testing)
#### proportions per cluster
testing.prop.sexes <- ggplot(data = testing.sex.distb.df.long, mapping = aes(x = CellType, y = count, fill = Sex)) +
  geom_bar(position="fill", stat="identity") +
  labs(x = "Cell Type", y = "Proportion") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 15),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 15),
        panel.grid = element_blank()) +
  NoLegend()

Rn7_VTA.testing@active.ident <- Rn7_VTA.testing$CellType # reset active ident to cell type

### Y gene counts
#### feature plot
testing.dim.plot.y <- FeaturePlot(Rn7_VTA.testing,
                                  features = "y.counts",
                                  split.by = "Sex",
                                  cols = c("lightgrey","blue"),
                                  keep.scale = "all")
#### violin
testing.vln.plot.y <- VlnPlot(Rn7_VTA.testing,
                              features = "y.counts",
                              split.by = "Sex",
                              pt.size = 0) +
  labs(x = "Cell Type", y = "ChrY Gene Expression") +
  NoLegend()

### Xist
#### feature plot
testing.dim.plot.xist <- FeaturePlot(Rn7_VTA.testing,
                                     features = "Xist.counts",
                                     split.by = "Sex",
                                     cols = c("lightgrey","red"),
                                     keep.scale = "all") 

#### violin
testing.vln.plot.xist <- VlnPlot(Rn7_VTA.testing,
                                 features = "Xist.counts",
                                 split.by = "Sex",
                                 pt.size = 0) +
  labs(x = "Cell Type", y = "Xist Gene Expression") +
  NoLegend()

### DEGs
testing.vln.plot.degs <- VlnPlot(Rn7_VTA.testing,
                                 features = c("ENSRNOG00000060617",
                                              "ENSRNOG00000065796",
                                              "Kdm5d",
                                              "Eif2s3y",
                                              "Kdm6a",
                                              "Cntnap2",
                                              "Crip1"),
                                 split.by = "Sex",
                                 pt.size = 0,
                                 stack = T,
                                 flip = T) +
  labs(x = "Cell Type") +
  theme(legend.position = "top")


### gene counts
#### distribution
testing.gene.count.distb <- ggplot(data = Rn7_VTA.testing@meta.data, mapping = aes(x = nFeature_RNA, fill = Neuronal)) +
  geom_density(alpha = 0.7) +
  labs(x = "Number of Genes", y = "Density") +
  theme_bw() +
  theme(text = element_text(size = 15),
        axis.text = element_text(size = 10),
        legend.position = c(0.9, 0.9),
        legend.title= element_blank()) +
  NoGrid()

# save outputs ####
### UMAPs
#### cell types
ggsave(plot = all.dimplot.clusters,
       file = "/path/to/sex_prediction_model/output/object_summary_plots/all.UMAP.celltypes.pdf",
       device = "pdf",
       height = 10,
       width = 12)
ggsave(plot = training.dimplot.clusters,
       file = "/path/to/sex_prediction_model/output/object_summary_plots/training.UMAP.celltypes.pdf",
       device = "pdf",
       height = 10,
       width = 12)
ggsave(plot = testing.dimplot.clusters,
       file = "/path/to/sex_prediction_model/output/object_summary_plots/testing.UMAP.celltypes.pdf",
       device = "pdf",
       height = 10,
       width = 12)
#### sexes
ggsave(plot = all.dimplot.sexes,
       file = "/path/to/sex_prediction_model/output/object_summary_plots/all.UMAP.sexes.pdf",
       device = "pdf",
       height = 10,
       width = 12)
ggsave(plot = training.dimplot.sexes,
       file = "/path/to/sex_prediction_model/output/object_summary_plots/training.UMAP.sexes.pdf",
       device = "pdf",
       height = 10,
       width = 12)
ggsave(plot = testing.dimplot.sexes,
       file = "/path/to/sex_prediction_model/output/object_summary_plots/testing.UMAP.sexes.pdf",
       device = "pdf",
       height = 10,
       width = 12)
#### ChrY expression
ggsave(plot = all.dim.plot.y,
       file = "/path/to/sex_prediction_model/output/object_summary_plots/all.UMAP.chrY.pdf",
       device = "pdf",
       height = 10,
       width = 20)
ggsave(plot = training.dim.plot.y,
       file = "/path/to/sex_prediction_model/output/object_summary_plots/training.UMAP.chrY.pdf",
       device = "pdf",
       height = 10,
       width = 20)
ggsave(plot = testing.dim.plot.y,
       file = "/path/to/sex_prediction_model/output/object_summary_plots/testing.UMAP.chrY.pdf",
       device = "pdf",
       height = 10,
       width = 20)
#### Xist expression
ggsave(plot = all.dim.plot.xist,
       file = "/path/to/sex_prediction_model/output/object_summary_plots/all.UMAP.xist.pdf",
       device = "pdf",
       height = 10,
       width = 20)
ggsave(plot = training.dim.plot.xist,
       file = "/path/to/sex_prediction_model/output/object_summary_plots/training.UMAP.xist.pdf",
       device = "pdf",
       height = 10,
       width = 20)
ggsave(plot = testing.dim.plot.xist,
       file = "/path/to/sex_prediction_model/output/object_summary_plots/testing.UMAP.xist.pdf",
       device = "pdf",
       height = 10,
       width = 20)

### Violin plots
#### ChrY expression
ggsave(plot = all.vln.plot.y,
       file = "/path/to/sex_prediction_model/output/object_summary_plots/all.vln.chrY.pdf",
       device = "pdf",
       height = 10,
       width = 12)
ggsave(plot = training.vln.plot.y,
       file = "/path/to/sex_prediction_model/output/object_summary_plots/training.vln.chrY.pdf",
       device = "pdf",
       height = 10,
       width = 12)
ggsave(plot = testing.vln.plot.y,
       file = "/path/to/sex_prediction_model/output/object_summary_plots/testing.vln.chrY.pdf",
       device = "pdf",
       height = 10,
       width = 12)
#### Xist expression
ggsave(plot = all.vln.plot.xist,
       file = "/path/to/sex_prediction_model/output/object_summary_plots/all.vln.xist.pdf",
       device = "pdf",
       height = 10,
       width = 12)
ggsave(plot = training.vln.plot.xist,
       file = "/path/to/sex_prediction_model/output/object_summary_plots/training.vln.xist.pdf",
       device = "pdf",
       height = 10,
       width = 12)
ggsave(plot = testing.vln.plot.xist,
       file = "/path/to/sex_prediction_model/output/object_summary_plots/testing.vln.xist.pdf",
       device = "pdf",
       height = 10,
       width = 12)
#### DEGs
ggsave(plot = all.vln.plot.degs,
       file = "/path/to/sex_prediction_model/output/object_summary_plots/all.vln.degs.pdf",
       device = "pdf",
       height = 10,
       width = 10)
ggsave(plot = training.vln.plot.degs,
       file = "/path/to/sex_prediction_model/output/object_summary_plots/training.vln.degs.pdf",
       device = "pdf",
       height = 10,
       width = 10)
ggsave(plot = testing.vln.plot.degs,
       file = "/path/to/sex_prediction_model/output/object_summary_plots/testing.vln.degs.pdf",
       device = "pdf",
       height = 10,
       width = 10)

### bar graphs
#### sex proportions
ggsave(plot = all.prop.sexes,
       file = "/path/to/sex_prediction_model/output/object_summary_plots/all.bar.sexes.pdf",
       device = "pdf",
       height = 10,
       width = 10)
ggsave(plot = training.prop.sexes,
       file = "/path/to/sex_prediction_model/output/object_summary_plots/training.bar.sexes.pdf",
       device = "pdf",
       height = 10,
       width = 10)
ggsave(plot = testing.prop.sexes,
       file = "/path/to/sex_prediction_model/output/object_summary_plots/testing.bar.sexes.pdf",
       device = "pdf",
       height = 10,
       width = 10)
#### ncells
ggsave(plot = all.ncells.clusters,
       file = "/path/to/sex_prediction_model/output/object_summary_plots/all.bar.ncells.pdf",
       device = "pdf",
       height = 10,
       width = 10)
ggsave(plot = training.ncells.clusters,
       file = "/path/to/sex_prediction_model/output/object_summary_plots/training.bar.ncells.pdf",
       device = "pdf",
       height = 10,
       width = 10)
ggsave(plot = testing.ncells.clusters,
       file = "/path/to/sex_prediction_model/output/object_summary_plots/testing.bar.ncells.pdf",
       device = "pdf",
       height = 10,
       width = 10)

### gene count distributions
ggsave(plot = all.gene.count.distb,
       file = "/path/to/sex_prediction_model/output/object_summary_plots/all.distb.genes.pdf",
       device = "pdf",
       height = 10,
       width = 10)
ggsave(plot = training.gene.count.distb,
       file = "/path/to/sex_prediction_model/output/object_summary_plots/training.distb.genes.pdf",
       device = "pdf",
       height = 10,
       width = 10)
ggsave(plot = testing.gene.count.distb,
       file = "/path/to/sex_prediction_model/output/object_summary_plots/testing.distb.genes.pdf",
       device = "pdf",
       height = 10,
       width = 10)


# sessionInfo ####
sessionInfo()