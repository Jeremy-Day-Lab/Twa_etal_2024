# Setup ####
# libraries
library(Seurat)
library(Boruta)
library(dplyr)
library(tidyr)
library(ggplot2)
library(plotly)
library(forcats)

# seed
set.seed(1234)

# data
### Training object
training.obj <- readRDS("/scratch/gtwa/Day/sex_prediction_model/Rn7_Xist_add/processed_data/2_create_splits/Rn7_VTA_training.RDS") 

## DEGs
### full results
degs.full <- read.csv("/scratch/gtwa/Day/sex_prediction_model/Rn7_Xist_add/processed_data/3_feature_selection_DEGs/DEGs_sex_cluster_VTA_full.csv") %>% 
  mutate(Cluster= factor(Cluster, levels = c("DA-Neuron",
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
                                             "Endothelial")))
### significant results
degs.sig <- degs.full %>% filter(p_val_adj < 0.05)
### significant results gene names
degs.sig.genes <- degs.sig %>% 
  pull(GeneName) %>% 
  unique()
### names of clusters
degs.sig.cluster <- degs.sig %>% 
  pull(Cluster) %>% 
  unique()

## Boruta
boruta.obj <- readRDS("/scratch/gtwa/Day/sex_prediction_model/Rn7_Xist_add/processed_data/4_feature_selection_Boruta/VTA_Boruta_max3000.RDS")
boruta.obj.fix <- TentativeRoughFix(boruta.obj)
selected.genes <- read.csv("/scratch/gtwa/Day/sex_prediction_model/Rn7_Xist_add/processed_data/4_feature_selection_Boruta/Boruta_final_decisions.csv") %>% 
  filter(finalDecision == "Confirmed") %>% 
  pull("variable")


# DEG plotting ####
Idents(training.obj) <- training.obj$CellType_Sex

## fill in lfc values for significant genes, across clusters
### initialize results list
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
### flatten list
lfc.df <- do.call(what = rbind,lfc.df)
### Create an ID column that is the rownames
lfc.df$ID <- rownames(lfc.df)
### Add gene names and cluster
lfc.df <- separate(lfc.df, col = ID, into = c("Cluster", "GeneName"), sep = "\\.", remove = FALSE, extra = "merge")
### filter for significant results
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
  mutate(Cluster= factor(Cluster, levels = c("DA-Neuron",
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
                                             "Endothelial") %>% rev()))

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

## save outs
ggsave(plot = degs.plot.hm,
       filename = "/scratch/gtwa/Day/sex_prediction_model/Rn7_Xist_add/plots/feature_selection_plotting/deg.sig.heatmap.pdf",
       height = 5,
       width = 10,
       units = "in",
       device = "pdf")
ggsave(plot = degs.plot.hm,
       filename = "/scratch/gtwa/Day/sex_prediction_model/Rn7_Xist_add/plots/feature_selection_plotting/deg.sig.heatmap.png",
       height = 5,
       width = 10,
       units = "in",
       dpi = 300,
       device = "png")
write.csv(lfc.df,
          file = "/scratch/gtwa/Day/sex_prediction_model/Rn7_Xist_add/processed_data/feature_selection_plotting/unfiltered_deg_results.csv",
          row.names = F)
write.csv(deg.df.full.coalesce,
          file = "/scratch/gtwa/Day/sex_prediction_model/Rn7_Xist_add/processed_data/feature_selection_plotting/lfc_sig_allclusters.csv", 
          row.names = F)


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
levels(training.obj$CellType) <- c("DA-Neuron",
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
Idents(training.obj) <- training.obj$CellType

top_10_summary <- df.long %>%
  group_by(name) %>%
  summarize(median_value = median(value, na.rm = TRUE)) %>%
  arrange(desc(median_value)) %>%
  slice_head(n = 10)

bottom_10_summary <- df.long %>%
  group_by(name) %>%
  summarize(median_value = median(value, na.rm = TRUE)) %>%
  arrange(desc(median_value)) %>%
  slice_tail(n = 10)


genes.vln <- VlnPlot(training.obj,
                     features = c("Xist",               # Xist most important female biased gene
                                  "ENSRNOG00000065796", # Xist proxy, second important female biased gene
                                  "ENSRNOG00000060617", # Uty, most important male biased gene
                                  "Ddx3",               # second most important male biased gene
                                  "Slc27a2",            # two least important genes
                                  "Trabd2b"),
                     split.by = "Sex",
                     cols = c("#F8766D", "#00BFC4"),
                     flip = T, stack = T) + theme(legend.position = "top")

## save outs
ggsave(plot = boruta.plot,
       filename = "/scratch/gtwa/Day/sex_prediction_model/Rn7_Xist_add/plots/feature_selection_plotting/boruta.decisions.pdf",
       height = 5,
       width = 10,
       units = "in",
       device = "pdf")
ggsave(plot = boruta.plot,
       filename = "/scratch/gtwa/Day/sex_prediction_model/Rn7_Xist_add/plots/feature_selection_plotting/boruta.decisions.png",
       height = 5,
       width = 10,
       units = "in",
       dpi = 300,
       device = "png")
ggsave(plot = genes.vln,
       filename = "/scratch/gtwa/Day/sex_prediction_model/Rn7_Xist_add/plots/feature_selection_plotting/example_gene_vln.pdf",
       height = 11,
       width = 7,
       units = "in",
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
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
# [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] forcats_1.0.0      plotly_4.10.4      ggplot2_3.5.1      tidyr_1.3.1        dplyr_1.1.4        Boruta_8.0.0       Seurat_5.2.1       SeuratObject_5.0.2 sp_2.1-4          
# 
# loaded via a namespace (and not attached):
#   [1] Rtsne_0.17             colorspace_2.1-0       deldir_2.0-4           ggridges_0.5.6         RcppHNSW_0.6.0         rstudioapi_0.16.0      spatstat.data_3.1-2   
# [8] farver_2.1.2           listenv_0.9.1          ggrepel_0.9.5          RSpectra_0.16-1        codetools_0.2-20       splines_4.2.0          polyclip_1.10-6       
# [15] spam_2.11-1            jsonlite_1.8.8         ica_1.0-3              cluster_2.1.6          png_0.1-8              uwot_0.2.2             shiny_1.8.1.1         
# [22] sctransform_0.4.1      spatstat.sparse_3.1-0  compiler_4.2.0         httr_1.4.7             Matrix_1.6-5           fastmap_1.2.0          lazyeval_0.2.2        
# [29] limma_3.54.2           cli_3.6.4              later_1.3.2            htmltools_0.5.8.1      tools_4.2.0            igraph_2.0.3           dotCall64_1.1-1       
# [36] gtable_0.3.5           glue_1.8.0             RANN_2.6.1             reshape2_1.4.4         Rcpp_1.0.14            scattermore_1.2        spatstat.univar_3.0-0 
# [43] vctrs_0.6.5            presto_1.0.0           spatstat.explore_3.3-1 nlme_3.1-165           progressr_0.14.0       lmtest_0.9-40          spatstat.random_3.3-1 
# [50] stringr_1.5.1          globals_0.16.3         mime_0.12              miniUI_0.1.1.1         lifecycle_1.0.4        irlba_2.3.5.1          goftest_1.2-3         
# [57] future_1.33.2          MASS_7.3-60.0.1        zoo_1.8-12             scales_1.3.0           ragg_1.3.2             promises_1.3.0         spatstat.utils_3.0-5  
# [64] parallel_4.2.0         RColorBrewer_1.1-3     reticulate_1.38.0      pbapply_1.7-2          gridExtra_2.3          stringi_1.8.4          fastDummies_1.7.3     
# [71] systemfonts_1.1.0      rlang_1.1.5            pkgconfig_2.0.3        matrixStats_1.3.0      lattice_0.22-6         ROCR_1.0-11            purrr_1.0.4           
# [78] tensor_1.5             labeling_0.4.3         patchwork_1.2.0        htmlwidgets_1.6.4      cowplot_1.1.3          tidyselect_1.2.1       parallelly_1.37.1     
# [85] RcppAnnoy_0.0.22       plyr_1.8.9             magrittr_2.0.3         R6_2.6.1               generics_0.1.3         pillar_1.10.1          withr_3.0.2           
# [92] fitdistrplus_1.2-1     survival_3.7-0         abind_1.4-5            tibble_3.2.1           future.apply_1.11.2    crayon_1.5.3           utf8_1.2.4            
# [99] KernSmooth_2.23-24     spatstat.geom_3.3-2    grid_4.2.0             data.table_1.17.0      digest_0.6.36          xtable_1.8-4           httpuv_1.6.15         
# [106] textshaping_0.4.0      munsell_0.5.1          viridisLite_0.4.2  
