# Setup ####
# Libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(tibble)
library(tidyr)

# Data
## Original NAc object
NAc_noXist <- readRDS("/scratch/gtwa/Day/sex_prediction_model/Rn7_Xist_add/raw_data/NAc_Combo_Integrated.RDS")
NAc_Xist <- readRDS("/scratch/gtwa/Day/sex_prediction_model/Rn7_Xist_add/processed_data/13_evaluation_NAc/NAc_Xist.RDS")

## NAc_Xist subsets
NAc_obj <- subset(NAc_noXist, cells = Cells(NAc_Xist))

# make metadata ordered factors
celltypes <- levels(NAc_obj$Combo_CellType) %>% rev()

NAc_obj$Combo_CellType <- factor(NAc_obj$Combo_CellType, levels = celltypes)
NAc_obj$Sex <- factor(NAc_obj$Sex, levels = c("Female", "Male"))

# Cluster summary plotting ####
# by cluster
NAc_obj@active.ident <- NAc_obj$Combo_CellType
all.dimplot.clusters <- DimPlot(NAc_obj, label = TRUE) + NoLegend()

# by sex
## all cells
NAc_obj@active.ident <- NAc_obj$Sex
all.dimplot.sex <- DimPlot(NAc_obj, shuffle = TRUE)

# Count cells
all.ncells.df <- table(NAc_obj$Combo_CellType) %>% # make a table of cell type counts
  data.frame() %>% # convert to dataframe
  rename("CellType" = "Var1",
         "Cells" = "Freq") %>% # rename columns
  mutate(CellType = factor(CellType, levels = celltypes)) # make cell type a factor

# CellType Cells
# 1         Drd1-MSN-1  5083
# 2         Drd1-MSN-2  2219
# 3         Drd2-MSN-1  5338
# 4         Drd2-MSN-2   490
# 5           Drd3-MSN   564
# 6           Grm8-MSN  2975
# 7          GABAergic  3079
# 8   Chat-Interneuron    84
# 9  Pvalb-Interneuron   724
# 10   Sst-Interneuron   403
# 11     Glutamatergic   680
# 12         Astrocyte  2893
# 13         Microglia  1503
# 14             Mural   128
# 15            Olig-1 10995
# 16    Polydendrocyte  2094

## save outputs
ggsave(plot = all.dimplot.clusters,
       file = "/scratch/gtwa/Day/sex_prediction_model/Rn7_Xist_add/plots/dataset_plotting/NAc_all.UMAP.celltypes.pdf",
       device = "pdf",
       height = 10,
       width = 10)

ggsave(plot = all.dimplot.sex,
       file = "/scratch/gtwa/Day/sex_prediction_model/Rn7_Xist_add/plots/dataset_plotting/NAc_all.UMAP.sex.pdf",
       device = "pdf",
       height = 10,
       width = 10)

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
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
# [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                 
# [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] tidyr_1.3.1        tibble_3.2.1       ggplot2_3.5.1      dplyr_1.1.4        Seurat_5.2.1       SeuratObject_5.0.2
# [7] sp_2.1-4          
# 
# loaded via a namespace (and not attached):
#   [1] Rtsne_0.17             colorspace_2.1-0       deldir_2.0-4           ggridges_0.5.6         RcppHNSW_0.6.0        
# [6] rstudioapi_0.16.0      spatstat.data_3.1-2    listenv_0.9.1          farver_2.1.2           ggrepel_0.9.5         
# [11] RSpectra_0.16-1        codetools_0.2-20       splines_4.2.0          polyclip_1.10-6        spam_2.11-1           
# [16] jsonlite_1.8.8         ica_1.0-3              cluster_2.1.6          png_0.1-8              uwot_0.2.2            
# [21] shiny_1.8.1.1          sctransform_0.4.1      spatstat.sparse_3.1-0  compiler_4.2.0         httr_1.4.7            
# [26] Matrix_1.6-5           fastmap_1.2.0          lazyeval_0.2.2         cli_3.6.4              later_1.3.2           
# [31] htmltools_0.5.8.1      tools_4.2.0            igraph_2.0.3           dotCall64_1.1-1        gtable_0.3.5          
# [36] glue_1.8.0             RANN_2.6.1             reshape2_1.4.4         Rcpp_1.0.14            scattermore_1.2       
# [41] spatstat.univar_3.0-0  vctrs_0.6.5            spatstat.explore_3.3-1 nlme_3.1-165           progressr_0.14.0      
# [46] lmtest_0.9-40          spatstat.random_3.3-1  stringr_1.5.1          globals_0.16.3         mime_0.12             
# [51] miniUI_0.1.1.1         lifecycle_1.0.4        irlba_2.3.5.1          goftest_1.2-3          future_1.33.2         
# [56] MASS_7.3-60.0.1        zoo_1.8-12             scales_1.3.0           ragg_1.3.2             promises_1.3.0        
# [61] spatstat.utils_3.0-5   parallel_4.2.0         RColorBrewer_1.1-3     reticulate_1.38.0      pbapply_1.7-2         
# [66] gridExtra_2.3          stringi_1.8.4          fastDummies_1.7.3      rlang_1.1.5            pkgconfig_2.0.3       
# [71] systemfonts_1.1.0      matrixStats_1.3.0      lattice_0.22-6         ROCR_1.0-11            purrr_1.0.4           
# [76] tensor_1.5             patchwork_1.2.0        htmlwidgets_1.6.4      labeling_0.4.3         cowplot_1.1.3         
# [81] tidyselect_1.2.1       parallelly_1.37.1      RcppAnnoy_0.0.22       plyr_1.8.9             magrittr_2.0.3        
# [86] R6_2.6.1               generics_0.1.3         pillar_1.10.1          withr_3.0.2            fitdistrplus_1.2-1    
# [91] survival_3.7-0         abind_1.4-5            future.apply_1.11.2    KernSmooth_2.23-24     spatstat.geom_3.3-2   
# [96] plotly_4.10.4          grid_4.2.0             data.table_1.17.0      digest_0.6.36          xtable_1.8-4          
# [101] httpuv_1.6.15          textshaping_0.4.0      munsell_0.5.1          viridisLite_0.4.2  