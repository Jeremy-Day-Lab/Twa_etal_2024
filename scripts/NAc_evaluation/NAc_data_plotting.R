# Setup ####
## libraries
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)

## set seed
set.seed(1234)

## data
### NAc testing data
testing.NAc <- readRDS("/data/project/daylab/2019-JD-0040/MCN_Code/Objects/NAc_Combo_Integrated.RDS")
### set order for cell types (neuronal -> nonneuronal)
testing.NAc$Combo_CellType <- factor(testing.NAc$Combo_CellType, levels = c("Drd1-MSN-1", "Drd1-MSN-2",
                                                                            "Drd2-MSN-1", "Drd2-MSN-2",
                                                                            "Drd3-MSN",
                                                                            "Grm8-MSN",
                                                                            "GABAergic",
                                                                            "Chat-Interneuron","Pvalb-Interneuron","Sst-Interneuron",
                                                                            "Glutamatergic",
                                                                            "Astrocyte", "Microglia", "Olig-1","Polydendrocyte", "Mural"))
### make celltype the active ident
testing.NAc@active.ident <- testing.NAc$Combo_CellType
### make frequency table of sex per cell type
sex.distb.df.long <- table(testing.NAc$Combo_CellType, testing.NAc$Sex) %>% # make a table of cell counts per cell type by sex
  as.data.frame.matrix() %>% # convert to dataframe
  rownames_to_column("CellType") %>% # make cell type a column
  mutate(CellType = factor(CellType, levels = c("Drd1-MSN-1", "Drd1-MSN-2",
                                                         "Drd2-MSN-1", "Drd2-MSN-2",
                                                         "Drd3-MSN",
                                                         "Grm8-MSN",
                                                         "GABAergic",
                                                         "Chat-Interneuron","Pvalb-Interneuron","Sst-Interneuron",
                                                         "Glutamatergic",
                                                         "Astrocyte", "Microglia", "Olig-1","Polydendrocyte", "Mural"))) %>% # make cell type a factor
  pivot_longer(cols = Female:Male, names_to = "Sex", values_to = "count") # pivot longer for plotting

# Plots ####
## cell type
### UMAP
umap.celltype <- DimPlot(testing.NAc,
                         label = T) + NoLegend()
### Freq table
table(testing.NAc$Combo_CellType)

## sex
### make sex the active identifier
testing.NAc@active.ident <- factor(testing.NAc$Sex)
### UMAP
umap.sex <- DimPlot(testing.NAc,
                    shuffle = T) + NoLegend()
### stacked bar chart
prop.sexes <- ggplot(data = sex.distb.df.long, mapping = aes(x = CellType, y = count, fill = Sex)) +
  geom_bar(position="fill", stat="identity") +
  labs(x = "Cell Type", y = "Proportion") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 15),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 15),
        panel.grid = element_blank()) +
  NoLegend()

## RNA features btwn dataset
testing.NAc@active.ident <- factor(testing.NAc$Dataset)
### nFeature
nFeature.vln <- VlnPlot(testing.NAc,
                        "nFeature_RNA",
                        pt.size = 0) +
  theme_bw() +
  labs(title = "", x = "Dataset", y = "nFeature") +
  theme(panel.grid = element_blank(), text = element_text(size=15)) +
  NoLegend()

nCount.vln <- VlnPlot(testing.NAc,
                      "nCount_RNA",
                      pt.size = 0) +
  theme_bw() +
  labs(title = "", x = "Dataset", y = "nCount") +
  theme(panel.grid = element_blank(), text = element_text(size=15)) +
  NoLegend()


# Save output ####
## cell type
### UMAP
ggsave(plot = umap.celltype,
       file = "/scratch/gtwa/Day/sex_prediction_model/output/NAc_evaluation/NAc_data_plot/UMAP_celltypes.pdf",
       height = 5,
       width = 5,
       device = "pdf")

## sex
### UMAP
ggsave(plot = umap.sex,
       file = "/scratch/gtwa/Day/sex_prediction_model/output/NAc_evaluation/NAc_data_plot/UMAP_sex.pdf",
       height = 5,
       width = 5,
       device = "pdf")
### stacked bar
ggsave(plot = prop.sexes,
       file = "/scratch/gtwa/Day/sex_prediction_model/output/NAc_evaluation/NAc_data_plot/proportion_sex.pdf",
       height = 5,
       width = 5,
       device = "pdf")

## RNA features
ggsave(plot = nFeature.vln,
       file = "/scratch/gtwa/Day/sex_prediction_model/output/NAc_evaluation/NAc_data_plot/Vln_nFeature.pdf",
       height = 5,
       width = 5,
       device = "pdf")
### stacked bar
ggsave(plot = nCount.vln,
       file = "/scratch/gtwa/Day/sex_prediction_model/output/NAc_evaluation/NAc_data_plot/Vln_nCount.pdf",
       height = 5,
       width = 5,
       device = "pdf")

# sessionInfo ####
sessionInfo()