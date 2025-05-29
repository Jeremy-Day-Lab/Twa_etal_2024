# Setup ####
# Libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(tibble)
library(tidyr)

# Data
## Original VTA object
Rn7_VTA_og <- readRDS("/scratch/gtwa/Day/sex_prediction_model/Rn7_Xist_add/raw_data/Rn7_VTA.RDS")
## training cells
Rn7_VTA.training_cells <- readRDS("/scratch/gtwa/Day/sex_prediction_model/Rn7_Xist_add/processed_data/2_create_splits/traincell_vector.RDS")
## testing cells
Rn7_VTA.testing_cells <- readRDS("/scratch/gtwa/Day/sex_prediction_model/Rn7_Xist_add/processed_data/2_create_splits/testcell_vector.RDS")

## Rn7_Xist subsets
Rn7_VTA <- subset(Rn7_VTA_og, cells = c(Rn7_VTA.training_cells, Rn7_VTA.testing_cells))
Rn7_VTA_training <- subset(Rn7_VTA, cells = Rn7_VTA.training_cells)
Rn7_VTA_testing <- subset(Rn7_VTA, cells = Rn7_VTA.testing_cells)

# make metadata ordered factors
celltypes <- c("DA-Neuron",
               "Glut-Neuron-1",
               "Glut-Neuron-2",
               "Glut-Neuron-3",
               "GABA-Neuron-1",
               "GABA-Neuron-2",
               "GABA-Neuron-3",
               "Olig-1",
               "Olig-2",
               "Olig-3",
               "Astrocyte",
               "Microglia",
               "Polydendrocyte",
               "OPC-Olig-1",
               "Mural",
               "Endothelial")

Rn7_VTA$CellType <- factor(Rn7_VTA$CellType, levels = celltypes)
Rn7_VTA_training$CellType <- factor(Rn7_VTA_training$CellType, levels = celltypes)
Rn7_VTA_testing$CellType <- factor(Rn7_VTA_testing$CellType, levels = celltypes)

Rn7_VTA$Sex <- factor(Rn7_VTA$Sex, levels = c("Female", "Male"))
Rn7_VTA_training$Sex <- factor(Rn7_VTA_training$Sex, levels = c("Female", "Male"))
Rn7_VTA_testing$Sex <- factor(Rn7_VTA_testing$Sex, levels = c("Female", "Male"))

# Cluster summary plotting ####
# by cluster
## all cells
Rn7_VTA@active.ident <- Rn7_VTA$CellType
all.dimplot.clusters <- DimPlot(Rn7_VTA, label = TRUE) + NoLegend()

## training
Rn7_VTA_training@active.ident <- Rn7_VTA_training$CellType
training.dimplot.clusters <- DimPlot(Rn7_VTA_training, label = TRUE) + NoLegend()

## testing
Rn7_VTA_testing@active.ident <- Rn7_VTA_testing$CellType
testing.dimplot.clusters <- DimPlot(Rn7_VTA_testing, label = TRUE) + NoLegend()

## save outputs
ggsave(plot = all.dimplot.clusters,
       file = "/scratch/gtwa/Day/sex_prediction_model/Rn7_Xist_add/plots/dataset_plotting/all.UMAP.celltypes.pdf",
       device = "pdf",
       height = 10,
       width = 10)
ggsave(plot = training.dimplot.clusters,
       file = "/scratch/gtwa/Day/sex_prediction_model/Rn7_Xist_add/plots/dataset_plotting/training.UMAP.celltypes.pdf",
       device = "pdf",
       height = 10,
       width = 10)
ggsave(plot = testing.dimplot.clusters,
       file = "/scratch/gtwa/Day/sex_prediction_model/Rn7_Xist_add/plots/dataset_plotting/testing.UMAP.celltypes.pdf",
       device = "pdf",
       height = 10,
       width = 10)

# by sex
## all cells
Rn7_VTA@active.ident <- Rn7_VTA$Sex
all.dimplot.sex <- DimPlot(Rn7_VTA, shuffle = TRUE)

## training
Rn7_VTA_training@active.ident <- Rn7_VTA_training$Sex
training.dimplot.sex <- DimPlot(Rn7_VTA_training, shuffle = TRUE)

## testing
Rn7_VTA_testing@active.ident <- Rn7_VTA_testing$Sex
testing.dimplot.sex <- DimPlot(Rn7_VTA_testing, shuffle = TRUE)

## save outputs
ggsave(plot = all.dimplot.sex,
       file = "/scratch/gtwa/Day/sex_prediction_model/Rn7_Xist_add/plots/dataset_plotting/all.UMAP.sex.pdf",
       device = "pdf",
       height = 10,
       width = 10)
ggsave(plot = training.dimplot.sex,
       file = "/scratch/gtwa/Day/sex_prediction_model/Rn7_Xist_add/plots/dataset_plotting/training.UMAP.sex.pdf",
       device = "pdf",
       height = 10,
       width = 10)
ggsave(plot = testing.dimplot.sex,
       file = "/scratch/gtwa/Day/sex_prediction_model/Rn7_Xist_add/plots/dataset_plotting/testing.UMAP.sex.pdf",
       device = "pdf",
       height = 10,
       width = 10)


# count ncells per cluster
## all cells
all.ncells.df <- table(Rn7_VTA$CellType) %>% # make a table of cell type counts
  data.frame() %>% # convert to dataframe
  rename("CellType" = "Var1",
         "Cells" = "Freq") %>% # rename columns
  mutate(CellType = factor(CellType, levels = c("DA-Neuron","Glut-Neuron-1", "Glut-Neuron-2", "Glut-Neuron-3", "GABA-Neuron-1", "GABA-Neuron-2", "GABA-Neuron-3",
                                                "Olig-1", "Olig-2", "Olig-3", "Astrocyte","Microglia","Polydendrocyte","OPC-Olig-1", "Mural", "Endothelial"))) # make cell type a factor
# CellType Cells
# 1       DA-Neuron   406
# 2   Glut-Neuron-1  2999
# 3   Glut-Neuron-2   697
# 4   Glut-Neuron-3   136
# 5   GABA-Neuron-1   894
# 6   GABA-Neuron-2   506
# 7   GABA-Neuron-3   320
# 8          Olig-1  6658
# 9          Olig-2  3120
# 10         Olig-3   724
# 11      Astrocyte  2559
# 12      Microglia  1118
# 13 Polydendrocyte  1560
# 14     OPC-Olig-1   249
# 15          Mural   139
# 16    Endothelial    64

train.ncells.df <- table(Rn7_VTA_training$CellType) %>% # make a table of cell type counts
  data.frame() %>% # convert to dataframe
  rename("CellType" = "Var1",
         "Cells" = "Freq") %>% # rename columns
  mutate(CellType = factor(CellType, levels = c("DA-Neuron","Glut-Neuron-1", "Glut-Neuron-2", "Glut-Neuron-3", "GABA-Neuron-1", "GABA-Neuron-2", "GABA-Neuron-3",
                                                "Olig-1", "Olig-2", "Olig-3", "Astrocyte","Microglia","Polydendrocyte","OPC-Olig-1", "Mural", "Endothelial"))) # make cell type a factor
# CellType Cells
# 1       DA-Neuron   285
# 2   Glut-Neuron-1  2101
# 3   Glut-Neuron-2   489
# 4   Glut-Neuron-3    96
# 5   GABA-Neuron-1   627
# 6   GABA-Neuron-2   355
# 7   GABA-Neuron-3   225
# 8          Olig-1  4661
# 9          Olig-2  2185
# 10         Olig-3   508
# 11      Astrocyte  1793
# 12      Microglia   783
# 13 Polydendrocyte  1092
# 14     OPC-Olig-1   175
# 15          Mural    98
# 16    Endothelial    46
test.ncells.df <- table(Rn7_VTA_testing$CellType) %>% # make a table of cell type counts
  data.frame() %>% # convert to dataframe
  rename("CellType" = "Var1",
         "Cells" = "Freq") %>% # rename columns
  mutate(CellType = factor(CellType, levels = c("DA-Neuron","Glut-Neuron-1", "Glut-Neuron-2", "Glut-Neuron-3", "GABA-Neuron-1", "GABA-Neuron-2", "GABA-Neuron-3",
                                                "Olig-1", "Olig-2", "Olig-3", "Astrocyte","Microglia","Polydendrocyte","OPC-Olig-1", "Mural", "Endothelial"))) # make cell type a factor
# CellType Cells
# 1       DA-Neuron   121
# 2   Glut-Neuron-1   898
# 3   Glut-Neuron-2   208
# 4   Glut-Neuron-3    40
# 5   GABA-Neuron-1   267
# 6   GABA-Neuron-2   151
# 7   GABA-Neuron-3    95
# 8          Olig-1  1997
# 9          Olig-2   935
# 10         Olig-3   216
# 11      Astrocyte   766
# 12      Microglia   335
# 13 Polydendrocyte   468
# 14     OPC-Olig-1    74
# 15          Mural    41
# 16    Endothelial    18


# sex distribution per cluster
## proportions per cluster
all.sex.distb.df.long <- table(Rn7_VTA$CellType, Rn7_VTA$Sex) %>% # make a table of cell counts per cell type by sex
  as.data.frame.matrix() %>% # convert to dataframe
  rownames_to_column("CellType") %>% # make cell type a column
  mutate(CellType = factor(CellType, levels = c("DA-Neuron","Glut-Neuron-1", "Glut-Neuron-2", "Glut-Neuron-3", "GABA-Neuron-1", "GABA-Neuron-2", "GABA-Neuron-3",
                                                "Olig-1", "Olig-2", "Olig-3", "Astrocyte","Microglia","Polydendrocyte","OPC-Olig-1", "Mural", "Endothelial") %>% rev())) %>% # make cell type a factor
  pivot_longer(cols = Female:Male, names_to = "Sex", values_to = "count") # pivot longer for plotting

## plot 
all.prop.sexes <- ggplot(data = all.sex.distb.df.long, mapping = aes(x = CellType, y = count, fill = Sex)) +
  geom_bar(position="fill", stat="identity") +
  labs(x = "Cell Type", y = "Proportion") +
  theme_bw() +
  theme(text = element_text(size = 15),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 15),
        panel.grid = element_blank()) +
  coord_flip() +
  NoLegend()

## save plot
ggsave(plot = all.prop.sexes,
       filename = "/scratch/gtwa/Day/sex_prediction_model/Rn7_Xist_add/plots/dataset_plotting/all.prop.sex.pdf",
       height = 10,
       width = 10,
       device = "pdf")

# sessionInfo ####
sessionInfo()
# R version 4.2.0 (2022-04-22)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Red Hat Enterprise Linux
# 
# Matrix products: default
# BLAS/LAPACK: /data/rc/apps/rc/software/FlexiBLAS/3.2.1-GCC-12.2.0/lib64/libflexiblas.so.3.2
# 
# locale:
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8   
# [6] LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C            
# [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] tidyr_1.3.1        tibble_3.2.1       ggplot2_3.5.1      dplyr_1.1.4        Seurat_5.2.1       SeuratObject_5.0.2 sp_2.1-4          
# 
# loaded via a namespace (and not attached):
#   [1] Rtsne_0.17             colorspace_2.1-0       deldir_2.0-4           ggridges_0.5.6         RcppHNSW_0.6.0         rstudioapi_0.16.0     
# [7] spatstat.data_3.1-2    listenv_0.9.1          farver_2.1.2           ggrepel_0.9.5          RSpectra_0.16-1        codetools_0.2-20      
# [13] splines_4.2.0          polyclip_1.10-6        spam_2.11-1            jsonlite_1.8.8         ica_1.0-3              cluster_2.1.6         
# [19] png_0.1-8              uwot_0.2.2             shiny_1.8.1.1          sctransform_0.4.1      spatstat.sparse_3.1-0  compiler_4.2.0        
# [25] httr_1.4.7             Matrix_1.6-5           fastmap_1.2.0          lazyeval_0.2.2         cli_3.6.4              later_1.3.2           
# [31] htmltools_0.5.8.1      tools_4.2.0            igraph_2.0.3           dotCall64_1.1-1        gtable_0.3.5           glue_1.8.0            
# [37] RANN_2.6.1             reshape2_1.4.4         Rcpp_1.0.14            scattermore_1.2        spatstat.univar_3.0-0  vctrs_0.6.5           
# [43] spatstat.explore_3.3-1 nlme_3.1-165           progressr_0.14.0       lmtest_0.9-40          spatstat.random_3.3-1  stringr_1.5.1         
# [49] globals_0.16.3         mime_0.12              miniUI_0.1.1.1         lifecycle_1.0.4        irlba_2.3.5.1          goftest_1.2-3         
# [55] future_1.33.2          MASS_7.3-60.0.1        zoo_1.8-12             scales_1.3.0           ragg_1.3.2             promises_1.3.0        
# [61] spatstat.utils_3.0-5   parallel_4.2.0         RColorBrewer_1.1-3     reticulate_1.38.0      pbapply_1.7-2          gridExtra_2.3         
# [67] stringi_1.8.4          fastDummies_1.7.3      rlang_1.1.5            pkgconfig_2.0.3        systemfonts_1.1.0      matrixStats_1.3.0     
# [73] lattice_0.22-6         ROCR_1.0-11            purrr_1.0.4            tensor_1.5             patchwork_1.2.0        htmlwidgets_1.6.4     
# [79] labeling_0.4.3         cowplot_1.1.3          tidyselect_1.2.1       parallelly_1.37.1      RcppAnnoy_0.0.22       plyr_1.8.9            
# [85] magrittr_2.0.3         R6_2.6.1               generics_0.1.3         pillar_1.10.1          withr_3.0.2            fitdistrplus_1.2-1    
# [91] survival_3.7-0         abind_1.4-5            future.apply_1.11.2    KernSmooth_2.23-24     spatstat.geom_3.3-2    plotly_4.10.4         
# [97] grid_4.2.0             data.table_1.17.0      digest_0.6.36          xtable_1.8-4           httpuv_1.6.15          textshaping_0.4.0     
# [103] munsell_0.5.1          viridisLite_0.4.2 
