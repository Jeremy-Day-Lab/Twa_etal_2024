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
training.object <- readRDS("/scratch/gtwa/Day/sex_prediction_model/data/celltype_models/Rn7_VTA_training_nonneuronal.RDS")
### testing object data
testing.object <- readRDS("/scratch/gtwa/Day/sex_prediction_model/data/celltype_models/Rn7_VTA_testing_nonneuronal.RDS")

## functions
Xist.count <- function(srat.object){
  # load gene expression count data
  Xist.counts <- as.data.frame(t(as.matrix(GetAssayData(object = srat.object, slot = "data", assay = "RNA")))) %>%
    pull("ENSRNOG00000065796") # pull out Xist (named as ENSRNOG00000065796)
  
  return(Xist.counts)
}
Xist.class.threshold <- function(srat.object){
  class.roc <- roc(srat.object$Sex,srat.object$Xist.counts, levels = c("Male", "Female"))
  threshold <- coords(class.roc, x = "best")
  return(threshold$threshold)
}
Xist.count.classify <- function(srat.object, threshold){
  classifications <- case_when(srat.object$Xist.counts > threshold ~ "Female",
                               TRUE ~ "Male")
  return(classifications)
}


# Prepare classifier ####
## Sum Xist counts
training.object$Xist.counts <- Xist.count(training.object)

## Xist summary stats
### distributions
#### overall
overall.distb <- ggplot(data = training.object@meta.data, mapping = aes(x = Xist.counts)) +
  geom_histogram(bins = 20, color = "black", fill = "lightgrey") +
  labs(title = "Xist gene count distribution", x = "Gene count", y = "Cell count") +
  theme_bw()

#### split by sex
wrap.distb <- ggplot(data = training.object@meta.data, mapping = aes(x = Xist.counts)) +
  geom_histogram(bins = 20, color = "black", fill = "lightgrey") +
  labs(title = "Xist gene count distribution", x = "Gene count", y = "Cell count") +
  facet_wrap(~ Sex) +
  theme_bw()

#### save plots
ggsave(plot = overall.distb,
       filename = "/scratch/gtwa/Day/sex_prediction_model/output/celltype_models/Xist.distb.training.nonneuronal.all.pdf",
       height = 5,
       width = 5)

ggsave(plot = wrap.distb,
       filename = "/scratch/gtwa/Day/sex_prediction_model/output/celltype_models/Xist.distb.training.nonneuronal.split.pdf",
       height = 5,
       width = 10)


### summary statistics
#### overall
overall.stats <- data.frame(min = min(training.object$Xist.counts), 
                            Q1 = quantile(training.object$Xist.counts, 0.25),
                            median = median(training.object$Xist.counts),
                            mean = mean(training.object$Xist.counts),
                            Q3 = quantile(training.object$Xist.counts, 0.75),
                            max = max(training.object$Xist.counts))
rownames(overall.stats) <- c("Overall")

#### split by sex
summary.stats <- training.object@meta.data %>% 
  group_by(Sex) %>% 
  summarize(
    min = min(Xist.counts),
    Q1 = quantile(Xist.counts, 0.25),
    median = median(Xist.counts),
    mean = mean(Xist.counts),
    Q3 = quantile(Xist.counts, 0.75),
    max = max(Xist.counts)
  ) %>% 
  data.frame() %>% 
  column_to_rownames(var = "Sex")

#### bind 
stats.df <- bind_rows(summary.stats, overall.stats)
#### write out
write.csv(stats.df, file = "/scratch/gtwa/Day/sex_prediction_model/output/celltype_models/Xist.distb.training.nonneuronal.sumstats.csv")

## Set classification threshold
### ROC curve
simple.roc <- roc(training.object$Sex,training.object$Xist.counts,
                  levels = c("Female", "Male"),
                  direction = ">"
)

### save plot 
roc.plot <- 
  ggroc(data = simple.roc, legacy.axes = F) +
  geom_abline(slope = 1, intercept = 1, linetype = "dashed") +
  theme_bw()
ggsave(filename = "/scratch/gtwa/Day/sex_prediction_model/output/celltype_models/Xist.training.nonneuronal.roc.pdf",
       plot = roc.plot,
       height = 5,
       width = 5)

### set threshold
threshold <- Xist.class.threshold(training.object)


# Make classification predictions ####
## sum chrY counts
testing.object$Xist.counts <- Xist.count(testing.object)
## make precition
testing.object$Sex.prediction <- Xist.count.classify(testing.object, 0)
testing.object$Sex.roc.prediction <- Xist.count.classify(testing.object, threshold)
## save predictions
testing.predictions <- testing.object@meta.data %>% select(Sex, Xist.counts, Sex.prediction, Sex.roc.prediction)
write.csv(testing.predictions, file = "/scratch/gtwa/Day/sex_prediction_model/output/celltype_models/Xist.class.predictions.nonneuronal.csv")


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