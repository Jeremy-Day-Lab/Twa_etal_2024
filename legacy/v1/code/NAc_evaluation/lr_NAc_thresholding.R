# Setup ####
## Libraries
library(SeuratObject)
library(Seurat)
library(caret)
library(ROCR)
library(pROC)
library(ggplot2)
library(dplyr)
library(tidyr)
library(magrittr)

## Data ####
### NAc testing data
testing.NAc <- readRDS("/data/project/daylab/2019-JD-0040/MCN_Code/Objects/NAc_Combo_Integrated.RDS")
### model genes
confirmed.genes <- read.csv("/scratch/gtwa/Day/sex_prediction_model/data/Boruta_final_decisions.csv") %>% 
  filter(finalDecision == "Confirmed") %>%
  pull(variable)
### models
lr_model <- readRDS("/scratch/gtwa/Day/sex_prediction_model/data/models/lr_model.RDS")

# NAc object preparation ####
## extract count matrices
testing.NAc.df <- as.data.frame(t(as.matrix(GetAssayData(object = testing.NAc, slot = "data",assay = "RNA")))) %>% 
  select(any_of(confirmed.genes))
### make sure cells are in the same order
all(names(Idents(testing.NAc)) == rownames(testing.NAc.df)) #TRUE 
### replace gene names that have "-" with "." to match gene name form the VTA models expect
colnames(testing.NAc.df) <- gsub("-", ".", colnames(testing.NAc.df))

## add sex identity columns
### character column
testing.NAc.df$Identity_bin_char <- factor(testing.NAc$Sex, levels = c("Male", "Female"))
### binary columns 
testing.NAc.df$Identity_bin <- ifelse(testing.NAc.df$Identity_bin_char  == "Female",
                                      1, #Females are 1
                                      0) #Males are 0
testing.NAc.df$Identity_bin <- as.factor(testing.NAc.df$Identity_bin)


# Make predictions ####
## logistic regression
lr_prediction <- predict(lr_model,testing.NAc.df[,1:285], type = "response")
lr_prediction_df <- data.frame(cell = names(lr_prediction),
                               prediction = lr_prediction,
                               binarized = ifelse(lr_prediction >0.5, 1, 0),
                               character = ifelse(lr_prediction >0.5, "Female", "Male") %>% as.factor(),
                               real = (as.numeric(testing.NAc.df$Identity_bin)-1),
                               celltype = testing.NAc$Combo_CellType,
                               GEM_Well = testing.NAc$GEM)
lr_prediction_df <- lr_prediction_df %>% 
  mutate(correct = binarized == real)

# Evaluate thresholding predictions ####
thresholding_stats_df <- data.frame(threshold = c(0,0.1, 0.25),
                                    ncell_included = c(nrow(lr_prediction_df),
                                                       lr_prediction_df %>% filter(prediction < 0.4 | prediction > 0.6) %>% nrow(),
                                                       lr_prediction_df %>% filter(prediction < 0.25 | prediction > 0.75) %>% nrow()),
                                    ncell_excluded = c(0,
                                                       nrow(lr_prediction_df) - lr_prediction_df %>% filter(prediction < 0.4 | prediction > 0.6) %>% nrow(),
                                                       nrow(lr_prediction_df) - lr_prediction_df %>% filter(prediction < 0.25 | prediction > 0.75) %>% nrow()),
                                    nfemale_included = c(lr_prediction_df %>% pull(real) %>% sum(),
                                                         lr_prediction_df %>% filter(prediction < 0.4 | prediction > 0.6) %>% pull(real) %>% sum(),
                                                         lr_prediction_df %>% filter(prediction < 0.25 | prediction > 0.75) %>% pull(real) %>% sum()),
                                    nfemale_excluded = c(0,
                                                         lr_prediction_df %>% filter(prediction > 0.4 & prediction < 0.6) %>% pull(real) %>% sum(),
                                                         lr_prediction_df %>% filter(prediction > 0.25 & prediction < 0.75) %>% pull(real) %>% sum()),
                                    acc_included = c(lr_prediction_df %>%
                                                       pull(correct) %>%
                                                       sum() / lr_prediction_df %>% nrow(),
                                                     lr_prediction_df %>%
                                                       filter(prediction < 0.4 | prediction > 0.6) %>%
                                                       pull(correct) %>% sum() / lr_prediction_df %>% filter(prediction < 0.4 | prediction > 0.6) %>% nrow(),
                                                     lr_prediction_df %>%
                                                       filter(prediction < 0.25 | prediction > 0.75) %>%
                                                       pull(correct) %>% sum() / lr_prediction_df %>% filter(prediction < 0.25 | prediction > 0.75) %>% nrow()),
                                    acc_excluded = c(NA,
                                                     lr_prediction_df %>%
                                                       filter(prediction > 0.4 & prediction < 0.6) %>%
                                                       pull(correct) %>% sum() / lr_prediction_df %>% filter(prediction > 0.4 & prediction < 0.6) %>% nrow(),
                                                     lr_prediction_df %>%
                                                       filter(prediction > 0.25 & prediction < 0.75) %>%
                                                       pull(correct) %>% sum() / lr_prediction_df %>% filter(prediction > 0.25 & prediction < 0.75) %>% nrow())
)


## Plot NAc data class prediction distribution
NAc.hist <- ggplot(data = lr_prediction_df, aes(x = prediction)) +
  geom_histogram(binwidth = 0.05, fill = "grey") +
  ylim(c(0,3800)) +
  labs(title = "NAc Data: Class Prediction Probabilities (bin size = 5)",
       x = "probability",
       y = "count") +
  theme_bw()

## Plot NAc data class predicitons with thresholds
NAc.hist.threshold <- ggplot(data = lr_prediction_df, aes(x = prediction)) +
  geom_histogram(binwidth = 0.05, fill = "grey") +
  ylim(c(0,3800)) +
  geom_vline(xintercept = 0.25, color = "red", linetype = "dotted") +
  geom_vline(xintercept = 0.75, color = "red", linetype = "dotted") +
  geom_vline(xintercept = 0.4, color = "blue", linetype = "dotted") +
  geom_vline(xintercept = 0.6, color = "blue", linetype = "dotted") +
  labs(title = "NAc Data: Class Prediction Probabilities (bin size = 5)",
       x = "probability",
       y = "count") +
  theme_bw()

### plot distribution of false classifications
NAc_false.hist.threshold <- lr_prediction_df %>% filter(correct == F) %>% 
  ggplot(aes(x = prediction)) +
  geom_histogram(binwidth = 0.05, fill = "grey") +
  ylim(c(0,3800)) +
  geom_vline(xintercept = 0.25, color = "red", linetype = "dotted") +
  geom_vline(xintercept = 0.75, color = "red", linetype = "dotted") +
  geom_vline(xintercept = 0.4, color = "blue", linetype = "dotted") +
  geom_vline(xintercept = 0.6, color = "blue", linetype = "dotted") +
  labs(title = "NAc Data: Class Prediction Probabilities (bin size = 5)",
       x = "probability",
       y = "count") +
  theme_bw()

NAc_false.hist.threshold.scaled <- lr_prediction_df %>% filter(correct == F) %>% 
  ggplot(aes(x = prediction)) +
  geom_histogram(binwidth = 0.05, fill = "grey") +
  #ylim(c(0,25000)) +
  geom_vline(xintercept = 0.25, color = "red", linetype = "dotted") +
  geom_vline(xintercept = 0.75, color = "red", linetype = "dotted") +
  geom_vline(xintercept = 0.4, color = "blue", linetype = "dotted") +
  geom_vline(xintercept = 0.6, color = "blue", linetype = "dotted") +
  labs(title = "NAc Data: Class Prediction Probabilities (bin size = 5)",
       x = "probability",
       y = "count") +
  theme_bw()

### plot distributions of true and false classifications
NAc_double.hist.threshold <- lr_prediction_df %>%
  ggplot(aes(x = prediction, fill = correct)) +
  geom_histogram(binwidth = 0.05, position = 'stack') +
  scale_fill_manual(values=c("black", "#00A651")) +
  geom_vline(xintercept = 0.25, color = "red", linetype = "dotted") +
  geom_vline(xintercept = 0.75, color = "red", linetype = "dotted") +
  geom_vline(xintercept = 0.4, color = "blue", linetype = "dotted") +
  geom_vline(xintercept = 0.6, color = "blue", linetype = "dotted") +
  labs(title = "NAc Data: Class Prediction Probabilities (bin size = 0.05)",
       x = "probability",
       y = "count") +
  theme_bw()

## Plot NAc data class binarized predicitons
NAc.hist.binarized <- ggplot(data = lr_prediction_df, aes(x = binarized)) +
  geom_histogram(binwidth = 0.05, fill = "grey") +
  ylim(c(0,25000)) +
  labs(title = "NAc Data: Class Prediction binarized (bin size = 5)",
       x = "probability",
       y = "count") +
  theme_bw()

NAc.hist.true <- ggplot(data = lr_prediction_df, aes(x = (as.numeric(real)-1))) +
  geom_histogram(binwidth = 0.05, fill = "grey") +
  ylim(c(0,25000)) +
  labs(title = "NAc Data: True Class Identities (bin size = 5)",
       x = "probability",
       y = "count") +
  theme_bw()

## Affected cell types
### NAc data
NAc_celltypes <- table(testing.NAc$CellType) %>% data.frame()
NAc_celltypes <- NAc_celltypes %>% 
  mutate(exclude.10 = lr_prediction_df %>% 
           filter(prediction > 0.4 & prediction < 0.6) %$% 
           table(.$celltype) %>% data.frame() %>% pull(Freq),
         exclude.25 = lr_prediction_df %>% 
           filter(prediction > 0.25 & prediction < 0.75) %$% 
           table(.$celltype) %>% data.frame() %>% pull(Freq))

NAc_celltype.hist <- lr_prediction_df %>% 
  mutate(neuronal = ifelse(celltype %in% c("Polydendrocyte", "Olig-1", "Mural", "Microglia", "Astrocyte"), "Non-Neuronal", "Neuronal")) %>% 
  ggplot(aes(x = prediction, fill = neuronal)) +
  geom_histogram(binwidth = 0.05, position = 'stack') +
  scale_fill_manual(values=c("#624185", "#FFA345")) +
  geom_vline(xintercept = 0.25, color = "red", linetype = "dotted") +
  geom_vline(xintercept = 0.75, color = "red", linetype = "dotted") +
  geom_vline(xintercept = 0.4, color = "blue", linetype = "dotted") +
  geom_vline(xintercept = 0.6, color = "blue", linetype = "dotted") +
  labs(title = "NAc Data: Class Prediction Probabilities (bin size = 0.05)",
       x = "probability",
       y = "count") +
  theme_bw()

# Save output ####
## accuracy w/ different thresholds
write.csv(thresholding_stats_df, file = "/scratch/gtwa/Day/sex_prediction_model/output/thresholding/LR_threshold_NAc_stats.csv")
## ncells per cell type excluded by decision buffer
write.csv(NAc_celltypes, file = "/scratch/gtwa/Day/sex_prediction_model/output/thresholding/LR_NAc_ncells.csv")

## NAc data predicted probability distribution w/o thresholds
ggsave(plot = NAc.hist,
       filename = "/scratch/gtwa/Day/sex_prediction_model/output/thresholding/LR_NAc_prediction_hist.pdf",
       device = "pdf",
       height = 6,
       width = 10)

## NAc data predicted probability distribution w/ thresholds
ggsave(plot = NAc.hist.threshold,
       filename = "/scratch/gtwa/Day/sex_prediction_model/output/thresholding/LR_NAc_prediction_threshold_hist.pdf",
       device = "pdf",
       height = 6,
       width = 10)

## NAc data false predictions distribution w/ threshold unscaled
ggsave(plot = NAc_false.hist.threshold,
       filename = "/scratch/gtwa/Day/sex_prediction_model/output/thresholding/LR_NAc_prediction_false_hist.pdf",
       device = "pdf",
       height = 6,
       width = 10)

## NAc data false predictions distribution w/ threshold scaled
ggsave(plot = NAc_false.hist.threshold.scaled,
       filename = "/scratch/gtwa/Day/sex_prediction_model/output/thresholding/LR_NAc_prediction_false_hist_scaled.pdf",
       device = "pdf",
       height = 6,
       width = 10)

## NAc data predicted probability distribution binarized
ggsave(plot = NAc.hist.binarized,
       filename = "/scratch/gtwa/Day/sex_prediction_model/output/thresholding/LR_NAc_prediction_binarized_hist.pdf",
       device = "pdf",
       height = 6,
       width = 10)

## NAc data predicted probability distribution true vals
ggsave(plot = NAc.hist.true,
       filename = "/scratch/gtwa/Day/sex_prediction_model/output/thresholding/LR_NAc_prediction_true_hist.pdf",
       device = "pdf",
       height = 6,
       width = 10)

## NAc data predicted probability distributions for true/false classifications
ggsave(plot = NAc_double.hist.threshold,
       filename = "/scratch/gtwa/Day/sex_prediction_model/output/thresholding/LR_NAc_prediction_true_false_hist.pdf",
       device = "pdf",
       height = 6,
       width = 10)

## NAc data predicted probability distributions for neuronal/non-neuronal classifications
ggsave(plot = NAc_celltype.hist,
       filename = "/scratch/gtwa/Day/sex_prediction_model/output/thresholding/LR_NAc_prediction_celltype_hist.pdf",
       device = "pdf",
       height = 6,
       width = 10)

# Session Info ####
sessionInfo()