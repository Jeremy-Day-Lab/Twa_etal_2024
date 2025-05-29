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
### testing data
Rn7_VTA_testing <- readRDS("/path/to/sex_prediction_model/data/Rn7_VTA_testing.RDS")
testing <- read.csv("/path/to/sex_prediction_model/data/test_count_data.csv", row.names = 1)
#### make character ID bin
testing$Identity_bin_char <- recode(testing$Identity_bin, `0` = "Male", `1` = "Female")
testing$Identity_bin_char <- factor(testing$Identity_bin_char, levels = c("Male", "Female"))
#### make ID bin a factor
testing$Identity_bin <- as.factor(testing$Identity_bin)
#### make cell type column
testing$CellType <- Rn7_VTA_testing@meta.data$CellType

### logistic regression model
lr_model <- readRDS("/path/to/sex_prediction_model/data/models/lr_model.RDS")

# Make test data set predictions ####
lr_prediction <- predict(lr_model,testing[,1:285], type = "response")
lr_prediction_df <- data.frame(cell = names(lr_prediction),
                               prediction = lr_prediction,
                               binarized = ifelse(lr_prediction >0.5, 1, 0),
                               character = ifelse(lr_prediction >0.5, "Female", "Male") %>% as.factor(),
                               real = (as.numeric(testing$Identity_bin)-1),
                               celltype = testing$CellType,
                               GEM_Well = Rn7_VTA_testing$GEM_Well)
lr_prediction_df <- lr_prediction_df %>% 
  mutate(correct = binarized == real)

# Evaluate prediction thresholds using test data ####
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

## Plot test data class prediction distribution
test.hist <- ggplot(data = lr_prediction_df, aes(x = prediction)) +
  geom_histogram(binwidth = 0.05, fill = "grey") +
  ylim(c(0,3800)) +
  labs(title = "Test Data: Class Prediction Probabilities (bin size = 5)",
       x = "probability",
       y = "count") +
  theme_bw()

## Plot test data class predicitons with thresholds
test.hist.threshold <- ggplot(data = lr_prediction_df, aes(x = prediction)) +
  geom_histogram(binwidth = 0.05, fill = "grey") +
  ylim(c(0,3800)) +
  geom_vline(xintercept = 0.25, color = "red", linetype = "dotted") +
  geom_vline(xintercept = 0.75, color = "red", linetype = "dotted") +
  geom_vline(xintercept = 0.4, color = "blue", linetype = "dotted") +
  geom_vline(xintercept = 0.6, color = "blue", linetype = "dotted") +
  labs(title = "Test Data: Class Prediction Probabilities (bin size = 5)",
       x = "probability",
       y = "count") +
  theme_bw()

### plot distribution of false classifications
test_false.hist.threshold <- lr_prediction_df %>% filter(correct == F) %>% 
  ggplot(aes(x = prediction)) +
  geom_histogram(binwidth = 0.05, fill = "grey") +
  ylim(c(0,3800)) +
  geom_vline(xintercept = 0.25, color = "red", linetype = "dotted") +
  geom_vline(xintercept = 0.75, color = "red", linetype = "dotted") +
  geom_vline(xintercept = 0.4, color = "blue", linetype = "dotted") +
  geom_vline(xintercept = 0.6, color = "blue", linetype = "dotted") +
  labs(title = "Test Data: Class Prediction Probabilities (bin size = 5)",
       x = "probability",
       y = "count") +
  theme_bw()

test_false.hist.threshold.scaled <- lr_prediction_df %>% filter(correct == F) %>% 
  ggplot(aes(x = prediction)) +
  geom_histogram(binwidth = 0.05, fill = "grey") +
  #ylim(c(0,3800)) +
  geom_vline(xintercept = 0.25, color = "red", linetype = "dotted") +
  geom_vline(xintercept = 0.75, color = "red", linetype = "dotted") +
  geom_vline(xintercept = 0.4, color = "blue", linetype = "dotted") +
  geom_vline(xintercept = 0.6, color = "blue", linetype = "dotted") +
  labs(title = "Test Data: Class Prediction Probabilities (bin size = 5)",
       x = "probability",
       y = "count") +
  theme_bw()

### plot distributions of true and false classifications
test_double.hist.threshold <- lr_prediction_df %>%
  ggplot(aes(x = prediction, fill = correct)) +
  geom_histogram(binwidth = 0.05, position = 'stack') +
  scale_fill_manual(values=c("black", "#00A651")) +
  geom_vline(xintercept = 0.25, color = "red", linetype = "dotted") +
  geom_vline(xintercept = 0.75, color = "red", linetype = "dotted") +
  geom_vline(xintercept = 0.4, color = "blue", linetype = "dotted") +
  geom_vline(xintercept = 0.6, color = "blue", linetype = "dotted") +
  labs(title = "Test Data: Class Prediction Probabilities (bin size = 0.05)",
       x = "probability",
       y = "count") +
  theme_bw()

## Plot test data class binarized predicitons
test.hist.binarized <- ggplot(data = lr_prediction_df, aes(x = binarized)) +
  geom_histogram(binwidth = 0.05, fill = "grey") +
  ylim(c(0,3800)) +
  labs(title = "Test Data: Class Prediction binarized (bin size = 5)",
       x = "probability",
       y = "count") +
  theme_bw()

test.hist.true <- ggplot(data = lr_prediction_df, aes(x = (as.numeric(real)-1))) +
  geom_histogram(binwidth = 0.05, fill = "grey") +
  ylim(c(0,3800)) +
  labs(title = "Test Data: True Class Identities (bin size = 5)",
       x = "probability",
       y = "count") +
  theme_bw()

## Affected cell types
### Test data
test_celltypes <- table(testing$CellType) %>% data.frame()
test_celltypes <- test_celltypes %>% 
  mutate(exclude.10 = lr_prediction_df %>% 
           filter(prediction > 0.4 & prediction < 0.6) %$% 
           table(.$celltype) %>% data.frame() %>% pull(Freq),
         exclude.25 = lr_prediction_df %>% 
           filter(prediction > 0.25 & prediction < 0.75) %$% 
           table(.$celltype) %>% data.frame() %>% pull(Freq))

test_celltype.hist <- lr_prediction_df %>% 
  mutate(neuronal = ifelse(celltype %in% c("Glut-Neuron-1", "GABA-Neuron-1", "Glut-Neuron-2", "GABA-Neuron-2", "DA-Neuron", "GABA-Neuron-3", "Glut-Neuron-3"), "Neuronal", "Non-neuronal")) %>% 
  ggplot(aes(x = prediction, fill = neuronal)) +
  geom_histogram(binwidth = 0.05, position = 'stack') +
  scale_fill_manual(values=c("#624185", "#FFA345")) +
  geom_vline(xintercept = 0.25, color = "red", linetype = "dotted") +
  geom_vline(xintercept = 0.75, color = "red", linetype = "dotted") +
  geom_vline(xintercept = 0.4, color = "blue", linetype = "dotted") +
  geom_vline(xintercept = 0.6, color = "blue", linetype = "dotted") +
  labs(title = "Test Data: Class Prediction Probabilities (bin size = 0.05)",
       x = "probability",
       y = "count") +
  theme_bw()

# Save output
## accuracy w/ different thresholds
write.csv(thresholding_stats_df, file = "/path/to/sex_prediction_model/output/thresholding/LR_threshold_test_stats.csv")
## ncells per cell type excluded by decision buffer
write.csv(test_celltypes, file = "/path/to/sex_prediction_model/output/thresholding/LR_test_ncells.csv")

## test data predicted probability distribution w/o thresholds
ggsave(plot = test.hist,
       filename = "/path/to/sex_prediction_model/output/thresholding/LR_test_prediction_hist.pdf",
       device = "pdf",
       height = 6,
       width = 10)

## test data predicted probability distribution w/ thresholds
ggsave(plot = test.hist.threshold,
       filename = "/path/to/sex_prediction_model/output/thresholding/LR_test_prediction_threshold_hist.pdf",
       device = "pdf",
       height = 6,
       width = 10)

## test data false predictions distribution w/ threshold unscaled
ggsave(plot = test_false.hist.threshold,
       filename = "/path/to/sex_prediction_model/output/thresholding/LR_test_prediction_false_hist.pdf",
       device = "pdf",
       height = 6,
       width = 10)

## test data false predictions distribution w/ threshold scaled
ggsave(plot = test_false.hist.threshold.scaled,
       filename = "/path/to/sex_prediction_model/output/thresholding/LR_test_prediction_false_hist_scaled.pdf",
       device = "pdf",
       height = 6,
       width = 10)

## test data predicted probability distribution binarized
ggsave(plot = test.hist.binarized,
       filename = "/path/to/sex_prediction_model/output/thresholding/LR_test_prediction_binarized_hist.pdf",
       device = "pdf",
       height = 6,
       width = 10)

## test data predicted probability distribution true vals
ggsave(plot = test.hist.true,
       filename = "/path/to/sex_prediction_model/output/thresholding/LR_test_prediction_true_hist.pdf",
       device = "pdf",
       height = 6,
       width = 10)

## NAc data predicted probability distributions for true/false classifications
ggsave(plot = test_double.hist.threshold,
       filename = "/path/to/sex_prediction_model/output/thresholding/LR_test_prediction_true_false_hist.pdf",
       device = "pdf",
       height = 6,
       width = 10)

## NAc data predicted probability distributions for neuronal/non-neuronal classifications
ggsave(plot = test_celltype.hist,
       filename = "/path/to/sex_prediction_model/output/thresholding/LR_test_prediction_celltype_hist.pdf",
       device = "pdf",
       height = 6,
       width = 10)

# Session info ####
sessionInfo()