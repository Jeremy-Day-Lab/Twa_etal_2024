# Setup ####
# libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyr)
library(infotheo)
library(mpmi)
library(forcats)

# set seed
set.seed(1234)

# data
## training object
Rn7_VTA.training <- readRDS("/scratch/gtwa/Day/sex_prediction_model/Rn7_Xist_add/processed_data/2_create_splits/Rn7_VTA_training.RDS")
Rn7_VTA.training@meta.data <- Rn7_VTA.training@meta.data %>%
  mutate(Neuronal = case_when(CellType %in% c("Glut-Neuron-1", "Glut-Neuron-2", "Glut-Neuron-3", "GABA-Neuron-1", "GABA-Neuron-2", "GABA-Neuron-3", "DA-Neuron") ~ "Neuronal",
                              CellType %in% c("Olig-1", "Olig-2", "Olig-3", "Astrocyte", "Polydendrocyte", "Microglia", "OPC-Olig-1", "Mural", "Endothelial") ~ "Non-neuronal"))

## genes
selected.genes <- read.csv("/scratch/gtwa/Day/sex_prediction_model/Rn7_Xist_add/processed_data/4_feature_selection_Boruta/Boruta_final_decisions.csv") %>% 
  filter(finalDecision == "Confirmed") %>% 
  pull("variable")

# mpmpi calculation ####
# extract count matrix
count.data <- as.data.frame(t(as.matrix(GetAssayData(object = Rn7_VTA.training,slot = "data",assay = "RNA")))) %>% 
  select(all_of(selected.genes))
# create cell lists for neuronal and non-neuronal cells
neuronal.cells <- Rn7_VTA.training@meta.data %>% filter(Neuronal == "Neuronal") %>% rownames()
nonneuronal.cells <- Rn7_VTA.training@meta.data %>% filter(Neuronal != "Neuronal") %>% rownames()

# create count dfs for neuronal and non-neuronal cells
count.data.neurons <- count.data[neuronal.cells,] 
count.data.nonneurons <- count.data[nonneuronal.cells,] 

# calculate mutual info
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

# merge to dataframe
minfo.training.df <- as.data.frame(do.call(rbind,minfo))
# modify dataframe
colnames(minfo.training.df) <- c("Gene","neuronal", "nonneuronal") #set column names to be informative
minfo.global.df <- minfo.training.df %>%
  mutate(neuronal = neuronal %>% as.numeric(),
         nonneuronal = nonneuronal %>% as.numeric())
# pivot longer
minfo.training.df.long <- minfo.training.df %>% 
  pivot_longer(cols = neuronal:nonneuronal,
               names_to = "cell.type",
               values_to = "mutual.info") %>% 
  mutate(cell.type = recode(cell.type, neuronal = "Neuronal", nonneuronal = "Non-neuronal")) %>% 
  mutate(cell.type = factor(cell.type, levels = c("Neuronal", "Non-neuronal")),
         mutual.info = mutual.info %>% as.numeric())

# write out
write.csv(minfo.training.df.long, "/scratch/gtwa/Day/sex_prediction_model/Rn7_Xist_add/processed_data/11_mutual_information/mutualinfo.training.df.csv")

# Plot mutual info distribution ####
# density
training.dens <- ggplot(data = minfo.training.df.long, aes(x = mutual.info, fill = cell.type)) +
  scale_fill_manual(values = c(Neuronal = "#624185", `Non-neuronal` = "#FFA345")) +
  geom_density(alpha = 0.7) +
  scale_x_continuous(trans = "log10") +
  labs(x = "Mutual Information", y = "Density", fill = "Cell Type") +
  theme_bw() +
  theme(text = element_text(size = 15), legend.position = c(0.2,0.9))

# save
ggsave(plot = training.dens,
       file = "/scratch/gtwa/Day/sex_prediction_model/Rn7_Xist_add/plots/11_mutual_information/training.dens.pdf",
       device = "pdf",
       height = 10,
       width = 10)

# cell type comparison ####
training.minfo.df.wide <- minfo.training.df.long %>% 
  pivot_wider(names_from = "cell.type",
              values_from = "mutual.info") %>% mutate(difference = Neuronal - `Non-neuronal`)

training.diff <- ggplot(training.minfo.df.wide, aes(x = fct_reorder(Gene, difference), y = difference)) +
  geom_point(shape = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Genes", y = "Mutual information (Neuronal - Non-neuronal)") +
  theme_bw() +
  theme(axis.text.x=element_blank(),
        panel.grid = element_blank())

# save
ggsave(plot = training.diff,
       file = "/scratch/gtwa/Day/sex_prediction_model/Rn7_Xist_add/plots/11_mutual_information/training.diff.rank.pdf",
       device = "pdf",
       height = 10,
       width = 10)
write.csv(training.minfo.df.wide, file = "/scratch/gtwa/Day/sex_prediction_model/Rn7_Xist_add/processed_data/11_mutual_information/training.minfo.wide.csv")


# Gene distributions ####
# sum number of expressed model genes in each cell type
nfeature.neurons <- data.frame(cell = row.names(count.data.neurons),
                                        cell.type = rep("Neuronal", nrow(count.data.neurons)),
                                        nFeature = rowSums(count.data.neurons != 0))
nfeature.nonneurons <- data.frame(cell = row.names(count.data.nonneurons),
                                           cell.type = rep("Non-neuronal", nrow(count.data.nonneurons)),
                                           nFeature = rowSums(count.data.nonneurons != 0))
nfeature <- rbind(nfeature.neurons,nfeature.nonneurons)


gene.distb <- ggplot(nfeature, aes(x = nFeature, fill = cell.type)) +
  geom_density(alpha = 0.7) +
  scale_fill_manual(values = c("Neuronal" = "#624185", "Non-neuronal" = "#ffa345")) +
  labs(x = "Number of genes", y = "Density") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.1,0.9),
        legend.title= element_blank(),
        text = element_text(size = 15))

# save output
ggsave(plot = gene.distb,
       filename = "/scratch/gtwa/Day/sex_prediction_model/Rn7_Xist_add/plots/11_mutual_information/training.model.nFeature.distb.pdf",
       height = 10,
       width = 10)

# sessionInfo ####
sessionInfo()