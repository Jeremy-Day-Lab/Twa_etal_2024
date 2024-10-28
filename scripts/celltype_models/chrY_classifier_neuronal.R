# Setup ####
## libraries
library(ggplot2)
library(dplyr)
library(tibble)
library(Seurat)
library(pROC)
library(caret)

## data
### training object data
training.object <- readRDS("/scratch/gtwa/Day/sex_prediction_model/data/celltype_models/Rn7_VTA_training_neuronal.RDS")
### testing object data
testing.object <- readRDS("/scratch/gtwa/Day/sex_prediction_model/data/celltype_models/Rn7_VTA_testing_neuronal.RDS")

### gene chromosome data
Rn7.gtf <- rtracklayer::import("/scratch/gtwa/Day/sex_prediction_model/data/Rattus_norvegicus.mRatBN7.2.105.gtf")
Rn7.gtf <- as.data.frame(Rn7.gtf)

## functions
y.count <- function(srat.object){
  # load gene expression count data
  srat.counts <- as.data.frame(t(as.matrix(GetAssayData(object = srat.object, slot = "data", assay = "RNA"))))
  
  # load gene annotations
  Rn7.gtf <- rtracklayer::import("/scratch/gtwa/Day/sex_prediction_model/data/Rattus_norvegicus.mRatBN7.2.105.gtf")
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
y.class.threshold <- function(srat.object){
  class.roc <- roc(srat.object$Sex,srat.object$y.counts)
  threshold <- coords(class.roc, x = "best")
  return(threshold$threshold)
}
y.count.classify <- function(srat.object, threshold){
  classifications <- case_when(srat.object$y.counts > threshold ~ "Male",
                               TRUE ~ "Female")
  return(classifications)
}


# Prepare classifiers ####
## Sum chrY counts
training.object$y.counts <- y.count(training.object)

## chrY summary stats
### distributions
#### overall
overall.distb <- ggplot(data = training.object@meta.data, mapping = aes(x = y.counts)) +
  geom_histogram(bins = 20, color = "black", fill = "lightgrey") +
  labs(title = "Chromosome Y gene count distribution", x = "Gene count", y = "Cell count") +
  theme_bw()

#### split by sex
wrap.distb <- ggplot(data = training.object@meta.data, mapping = aes(x = y.counts)) +
  geom_histogram(bins = 20, color = "black", fill = "lightgrey") +
  labs(title = "Chromosome Y gene count distribution", x = "Gene count", y = "Cell count") +
  facet_wrap(~ Sex) +
  theme_bw()

#### save plots
ggsave(plot = overall.distb,
       filename = "/scratch/gtwa/Day/sex_prediction_model/output/celltype_models/Ydistb.training.neuronal.all.pdf",
       height = 5,
       width = 5)

ggsave(plot = wrap.distb,
       filename = "/scratch/gtwa/Day/sex_prediction_model/output/celltype_models/Ydistb.training.neuronal.split.pdf",
       height = 5,
       width = 10)


### summary statistics
#### overall
overall.stats <- data.frame(min = min(training.object$y.counts), 
                            Q1 = quantile(training.object$y.counts, 0.25),
                            median = median(training.object$y.counts),
                            mean = mean(training.object$y.counts),
                            Q3 = quantile(training.object$y.counts, 0.75),
                            max = max(training.object$y.counts))
rownames(overall.stats) <- c("Overall")

#### split by sex
summary.stats <- training.object@meta.data %>% 
  group_by(Sex) %>% 
  summarize(
    min = min(y.counts),
    Q1 = quantile(y.counts, 0.25),
    median = median(y.counts),
    mean = mean(y.counts),
    Q3 = quantile(y.counts, 0.75),
    max = max(y.counts)
  ) %>% 
  data.frame() %>% 
  column_to_rownames(var = "Sex")

#### bind 
stats.df <- bind_rows(summary.stats, overall.stats)
#### write out
write.csv(stats.df, file = "/scratch/gtwa/Day/sex_prediction_model/output/celltype_models/Ydistb.training.neuronal.sumstats.csv")

## Set classification threshold
### ROC curve
simple.roc <- roc(training.object$Sex,training.object$y.counts)

### save plot 
roc.plot <- ggroc(data = simple.roc) +
  geom_abline(slope = 1, intercept = 1, linetype = "dashed") +
  theme_bw()
ggsave(filename = "/scratch/gtwa/Day/sex_prediction_model/output/celltype_models/Ydistb.training.neuronal.roc.pdf",
       plot = roc.plot,
       height = 5,
       width = 5)

### set threshold
threshold <- y.class.threshold(training.object)


# Make classification predictions ####
## sum chrY counts
testing.object$y.counts <- y.count(testing.object)
## make precition
testing.object$Sex.prediction <- y.count.classify(testing.object, 0)
testing.object$Sex.roc.prediction <- y.count.classify(testing.object, threshold)
## save predictions
testing.predictions <- testing.object@meta.data %>% select(Sex, y.counts, Sex.prediction, Sex.roc.prediction)
write.csv(testing.predictions, file = "/scratch/gtwa/Day/sex_prediction_model/output/celltype_models/Yclass.predictions.neuronal.csv")


# Evaluate performance ####
## make Sex and predicted Sex a factor
testing.object@meta.data <- testing.object@meta.data %>% mutate(Sex = factor(testing.object$Sex, levels = c("Male", "Female")),
                                                                Sex.prediction = factor(testing.object$Sex.prediction, levels = c("Male", "Female")),
                                                                Sex.roc.prediction = factor(testing.object$Sex.roc.prediction, levels = c("Male", "Female")))

print("Threshold 0 predictions")
confusionMatrix(testing.object$Sex.prediction, testing.object$Sex)

print(paste0("Threshold ", threshold, " predictions"))
confusionMatrix(testing.object$Sex.roc.prediction, testing.object$Sex)


# sessionInfo ####
sessionInfo()