# Setup ####
## libraries
library(Seurat)
library(Boruta)
library(dplyr)
library(tidyr)
library(ggplot2)
library(plotly)
library(forcats)


## seed
set.seed(1234)

## data
### Training object
training.obj <- readRDS("/path/to/sex_prediction_model/data/Rn7_VTA_training.RDS") 
### Testing object
testing.obj <- readRDS("/path/to/sex_prediction_model/data/Rn7_VTA_testing.RDS")

### DEGs
#### full results set
degs.full <- read.csv("/path/to/sex_prediction_model/data/DEGs_sex_cluster_VTA_full.csv") %>% 
  mutate(Cluster= factor(Cluster, levels = c("DA-Neuron","Glut-Neuron-1", "Glut-Neuron-2", "Glut-Neuron-3", "GABA-Neuron-1", "GABA-Neuron-2", "GABA-Neuron-3",
                                                "Olig-1", "Olig-2", "Olig-3", "Astrocyte","Microglia","Polydendrocyte","OPC-Olig-1", "Mural", "Endothelial")))
#### significant results set
degs.sig <- degs.full %>% filter(p_val_adj < 0.05)
#### significant gene names
degs.sig.genes <- degs.sig %>% 
  pull(GeneName) %>% 
  unique()
#### names of clusters
degs.sig.cluster <- degs.sig %>% 
  pull(Cluster) %>% 
  unique()
### Boruta objet
boruta.obj <- readRDS("/path/to/sex_prediction_model/data/VTA_Boruta_max2000.RDS")
boruta.obj.fix <- TentativeRoughFix(boruta.obj)
#### genes
selected.genes <- read.csv("/path/to/sex_prediction_model/data/Boruta_final_decisions.csv") %>% 
  filter(finalDecision == "Confirmed") %>% 
  pull("variable")



# DEGs ####
## fill in lfc values for significant genes, across clusters
### calculate lfc values (note: there is a known issue in Seurat 4.3.0.1 where `FoldChange()` does not return consistent results with `FindMarkers()`. `FindMarkers()` values are correct.)
Idents(training.obj) <- training.obj$CellType_Sex

#### initialize results list
lfc.df <- vector(mode = "list",length = 16)
names(lfc.df) <- levels(training.obj$CellType)
for(i in names(lfc.df)){
  lfc.df[[i]] <- FindMarkers(object = training.obj,
                               ident.1 = paste0(i,"_Female"),
                               ident.2 = paste0(i,"_Male"),
                               logfc.threshold = 0,
                               min.pct = 0,
                               min.cells.group = 0)
}
#### flatten list
lfc.df <- do.call(what = rbind,lfc.df)

#### Create an ID column that is the rownames
lfc.df$ID <- rownames(lfc.df)

#### Add gene names and cluster
lfc.df <- separate(lfc.df, col = ID, into = c("Cluster", "GeneName"), sep = "\\.", remove = FALSE, extra = "merge")
#### filter for significant results
lfc.df <- lfc.df %>% 
  filter(GeneName %in% degs.sig.genes)

### make complete df
deg.df.full <- expand.grid(Cluster = levels(degs.sig.cluster),
                           GeneName = degs.sig.genes)
deg.df.full <- left_join(x = deg.df.full,
                         y = degs.sig %>% select(Cluster, GeneName, avg_log2FC),
                         by = c("Cluster", "GeneName"),
                         keep = F,
                         relationship = "one-to-one")
deg.df.full <- left_join(deg.df.full,
                         y = lfc.df %>% select(Cluster, GeneName, avg_log2FC),
                         by = c("Cluster", "GeneName"),
                         keep = F,
                         relationship = "one-to-one")

deg.df.full.coalesce <- deg.df.full %>% 
  mutate(avg_log2FC = coalesce(avg_log2FC.x, avg_log2FC.y)) %>% 
  select(Cluster, GeneName, avg_log2FC) %>% 
  mutate(Cluster= factor(Cluster, levels = c("DA-Neuron","Glut-Neuron-1", "Glut-Neuron-2", "Glut-Neuron-3", "GABA-Neuron-1", "GABA-Neuron-2", "GABA-Neuron-3",
                                             "Olig-1", "Olig-2", "Olig-3", "Astrocyte","Microglia","Polydendrocyte","OPC-Olig-1", "Mural", "Endothelial") %>% rev()))

## plot heatmap
degs.plot.hm <- ggplot(deg.df.full.coalesce , aes(x = Cluster, y = fct_reorder(GeneName, avg_log2FC, mean), fill = avg_log2FC)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "#00BFC4",
    mid = "white",
    high = "#F8766D",
    midpoint = 0,
    limits = c(-2.5,2.5),
    oob = scales::squish) +
  labs(x = "Cell Type", y = "Gene", fill = "Log2 FC") +
  coord_flip() +
  theme(axis.text.x = element_blank(),
        panel.background = element_blank())
### interactive version
degs.plotly <- ggplotly(degs.plot.hm)

# Boruta ####
#ggplotting code from https://stackoverflow.com/questions/73415232/how-to-use-ggplot2-to-plot-box-plots-from-borutas-results-in-r
## object manipulation
### create list of importance histories for each gene
lz <- lapply(1:ncol(boruta.obj.fix$ImpHistory),
             function(i) boruta.obj.fix$ImpHistory[is.finite(boruta.obj.fix$ImpHistory[,i]),i])
### set names for each gene list
names(lz) <- colnames(boruta.obj.fix$ImpHistory) 
### set order for list items to be decreasing by median importance
ii <- order(sapply(lz,stats::median))
lz[ii] -> lz
### collapse to data frame
lz_df <- do.call(rbind.data.frame, lz)
df <- as.data.frame(t(lz_df))
names(df) <- names(lz)
rownames(df) <- NULL
### make vector of colors
generateCol<-function(x,colCode,col,numShadow){
  #Checking arguments
  if(is.null(col) & length(colCode)!=4)
    stop('colCode should have 4 elements.')
  #Generating col
  if(is.null(col)){
    rep(colCode[4],length(x$finalDecision)+numShadow)->cc
    cc[c(x$finalDecision=='Confirmed',rep(FALSE,numShadow))]<-colCode[1]
    cc[c(x$finalDecision=='Tentative',rep(FALSE,numShadow))]<-colCode[2]
    cc[c(x$finalDecision=='Rejected',rep(FALSE,numShadow))]<-colCode[3]
    col=cc
  }
  return(col)
}
boruta.col <- generateCol(boruta.obj.fix, c('confirmed','tentative','rejected','shadow'), NULL, 3)
boruta.col <- boruta.col[ii]
color.match <- data.frame(name = names(lz),
                          color = factor(boruta.col, levels = c("confirmed", "rejected", "tentative", "shadow")))

df.long <- df %>% 
  pivot_longer(everything()) %>% 
  mutate(name = factor(name, levels = names(df)))

df.long <- left_join(df.long,
                          color.match,
                          by = "name",
                          keep = F,
                          relationship = "many-to-one")

# gene.highlights <- c("shadowMin", "shadowMean", "shadowMax")

boruta.plot <- ggplot(df.long, aes(x = fct_reorder(name, value, median), y = value)) +
  geom_boxplot(mapping = aes(col = color), fill = "white",
               outlier.alpha = 0.5,
               outlier.size = 0.5,
               outlier.shape = NA) +
  # scale_x_discrete(breaks = gene.highlights,
  #                  labels = c("Shadow Min", "Shadow Mean", "Shadow Max")) +
  scale_color_manual(values = c("shadow" = "blue", "confirmed" = "green", "rejected" = "red", "tentative" = "yellow")) +
  labs(x = "Gene", y = "Importance", color = "Decision") +
  theme_bw() + 
  theme(legend.position = c(0.15,0.85),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_blank(),
        text = element_text(size = 15))

## plot example genes
Idents(training.obj) <- training.obj$CellType

genes.vln <- VlnPlot(training.obj,
                     features = c("ENSRNOG00000060617", # Uty, most important male biased gene
                                  "Ddx3",               # second most important male biased gene
                                  "ENSRNOG00000065796", # Xist, most important female biased gene
                                  "Kdm6a",              # second most important female biased gene
                                  "Fam210a",            # two least important genes
                                  "Trabd2b"),
                     split.by = "Sex",
                     cols = c("#F8766D", "#00BFC4"),
                     flip = T, stack = T) + theme(legend.position = "top")


# Gene distributions ####
# Pre-process data
## training
### make neuronal and non-neuronal metadata
training.obj@meta.data <- training.obj@meta.data %>%
  mutate(Neuronal = case_when(CellType %in% c("Glut-Neuron-1", "Glut-Neuron-2", "Glut-Neuron-3", "GABA-Neuron-1", "GABA-Neuron-2", "GABA-Neuron-3", "DA-Neuron") ~ "Neuronal",
                              CellType %in% c("Olig-1", "Olig-2", "Olig-3", "Astrocyte", "Polydendrocyte", "Microglia", "OPC-Olig-1", "Mural", "Endothelial") ~ "Non-neuronal"))
training.obj@meta.data <- training.obj@meta.data %>% 
  mutate(Neuronal = factor(Neuronal, levels = c("Neuronal", "Non-neuronal")))
### extract count matrix
training.count.data <- as.data.frame(t(as.matrix(GetAssayData(object = training.obj,slot = "data",assay = "RNA")))) %>% 
  select(all_of(selected.genes))
### create cell lists for neuronal and non-neuronal cells
training.neuronal.cells <- training.obj@meta.data %>% filter(Neuronal == "Neuronal") %>% rownames()
training.nonneuronal.cells <- training.obj@meta.data %>% filter(Neuronal != "Neuronal") %>% rownames()
### create count dfs for neuronal and non-neuronal cells
training.count.data.neurons <- training.count.data[training.neuronal.cells,] 
training.count.data.nonneurons <- training.count.data[training.nonneuronal.cells,]
### nFeatures
training.nfeature.neurons <- data.frame(cell = row.names(training.count.data.neurons),
                                        cell.type = rep("Neuronal", nrow(training.count.data.neurons)),
                                        nFeature = rowSums(training.count.data.neurons != 0))
training.nfeature.nonneurons <- data.frame(cell = row.names(training.count.data.nonneurons),
                                           cell.type = rep("Non-neuronal", nrow(training.count.data.nonneurons)),
                                           nFeature = rowSums(training.count.data.nonneurons != 0))
training.nfeature <- rbind(training.nfeature.neurons,training.nfeature.nonneurons)


training.gene.distb <- ggplot(training.nfeature, aes(x = nFeature, fill = cell.type)) +
  geom_density(alpha = 0.7) +
  scale_fill_manual(values = c("Neuronal" = "#624185", "Non-neuronal" = "#ffa345")) +
  labs(x = "Number of genes", y = "Density") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.1,0.9),
        legend.title= element_blank(),
        text = element_text(size = 15))

## testing
### make neuronal and non-neuronal metadata
testing.obj@meta.data <- testing.obj@meta.data %>%
  mutate(Neuronal = case_when(CellType %in% c("Glut-Neuron-1", "Glut-Neuron-2", "Glut-Neuron-3", "GABA-Neuron-1", "GABA-Neuron-2", "GABA-Neuron-3", "DA-Neuron") ~ "Neuronal",
                              CellType %in% c("Olig-1", "Olig-2", "Olig-3", "Astrocyte", "Polydendrocyte", "Microglia", "OPC-Olig-1", "Mural", "Endothelial") ~ "Non-neuronal"))
testing.obj@meta.data <- testing.obj@meta.data %>% 
  mutate(Neuronal = factor(Neuronal, levels = c("Neuronal", "Non-neuronal")))
### extract count matrix
testing.count.data <- as.data.frame(t(as.matrix(GetAssayData(object = testing.obj,slot = "data",assay = "RNA")))) %>% 
  select(all_of(selected.genes))
### create cell lists for neuronal and non-neuronal cells
testing.neuronal.cells <- testing.obj@meta.data %>% filter(Neuronal == "Neuronal") %>% rownames()
testing.nonneuronal.cells <- testing.obj@meta.data %>% filter(Neuronal != "Neuronal") %>% rownames()
### create count dfs for neuronal and non-neuronal cells
testing.count.data.neurons <- testing.count.data[testing.neuronal.cells,] 
testing.count.data.nonneurons <- testing.count.data[testing.nonneuronal.cells,] 
### nFeatures
testing.nfeature.neurons <- data.frame(cell = row.names(testing.count.data.neurons),
                                       cell.type = rep("Neuronal", nrow(testing.count.data.neurons)),
                                       nFeature = rowSums(testing.count.data.neurons != 0))
testing.nfeature.nonneurons <- data.frame(cell = row.names(testing.count.data.nonneurons),
                                          cell.type = rep("Non-neuronal", nrow(testing.count.data.nonneurons)),
                                          nFeature = rowSums(testing.count.data.nonneurons != 0))
testing.nfeature <- rbind(testing.nfeature.neurons,testing.nfeature.nonneurons)


testing.gene.distb <- ggplot(testing.nfeature, aes(x = nFeature, fill = cell.type)) +
  geom_density(alpha = 0.7) +
  scale_fill_manual(values = c("Neuronal" = "#624185", "Non-neuronal" = "#ffa345")) +
  labs(x = "Number of genes", y = "Density") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.1,0.9),
        legend.title= element_blank(),
        text = element_text(size = 15))

# Save outputs ####
## DEGs
ggsave(plot = degs.plot.hm,
       filename = "/path/to/sex_prediction_model/output/feature_selection/deg.sig.heatmap.pdf",
       height = 5,
       width = 10,
       units = "in",
       device = "pdf")
htmlwidgets::saveWidget(as_widget(degs.plotly), "/path/to/sex_prediction_model/output/feature_selection/deg.sig.heatmap.html")
write.csv(lfc.df,
          file = "/path/to/sex_prediction_model/output/feature_selection/unfiltered_deg_results.csv",
          row.names = F)
write.csv(deg.df.full.coalesce,
          file = "/path/to/sex_prediction_model/output/feature_selection/lfc_sig_allclusters.csv", 
          row.names = F)

## Boruta
ggsave(plot = boruta.plot,
       filename = "/path/to/sex_prediction_model/output/feature_selection/boruta.decisions.pdf",
       height = 5,
       width = 10,
       units = "in",
       device = "pdf")

ggsave(plot = genes.vln,
       filename = "/path/to/sex_prediction_model/output/feature_selection/example_gene_vln.pdf",
       height = 11,
       width = 7,
       units = "in",
       device = "pdf")

## model gene expression
ggsave(plot = training.gene.distb,
       filename = "/path/to/sex_prediction_model/output/feature_selection/training.model.nFeature.distb.pdf",
       height = 10,
       width = 10)

ggsave(plot = testing.gene.distb,
       filename = "/path/to/sex_prediction_model/output/feature_selection/testing.model.nFeature.distb.pdf",
       height = 10,
       width = 10)

# SessionInfo ####
sessionInfo()