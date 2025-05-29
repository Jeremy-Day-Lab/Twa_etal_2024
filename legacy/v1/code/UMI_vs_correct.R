# Setup ####
## libraries
library(Seurat)
library(caret)
library(randomForest)
library(kernlab)
library(RSNNS)
library(ROCR)
library(pROC)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

## set seed
set.seed(1234)

## functions
### Chr Y classification
y.count <- function(srat.object){
  # load gene expression count data
  srat.counts <- as.data.frame(t(as.matrix(GetAssayData(object = srat.object, slot = "data", assay = "RNA"))))
  
  # load gene annotations
  Rn7.gtf <- rtracklayer::import("/path/to/sex_prediction_model/data/Rattus_norvegicus.mRatBN7.2.105.gtf")
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
y.count.classify <- function(srat.object, threshold){
  srat.object$y.counts <- y.count(srat.object)
  classifications <- case_when(srat.object$y.counts > threshold ~ "Male",
                               TRUE ~ "Female")
  return(classifications)
}

### Xist classification
Xist.count <- function(srat.object){
  # load gene expression count data
  Xist.counts <- as.data.frame(t(as.matrix(GetAssayData(object = srat.object, slot = "data", assay = "RNA")))) %>%
    pull("ENSRNOG00000065796") # pull out Xist (named as ENSRNOG00000065796)
  
  return(Xist.counts)
}
Xist.count.classify <- function(srat.object, threshold){
  srat.object$Xist.counts <- Xist.count(srat.object)
  classifications <- case_when(srat.object$Xist.counts > threshold ~ "Female",
                               TRUE ~ "Male")
  return(classifications)
}

## data
### models
svm_model <- readRDS("/path/to/sex_prediction_model/data/models/svm_rad_model.RDS")
rf_model <- readRDS("/path/to/sex_prediction_model/data/models/rf_model.RDS")
mlp_model <- readRDS("/path/to/sex_prediction_model/data/models/mlpwdml_model.RDS")
lr_model <- readRDS("/path/to/sex_prediction_model/data/models/lr_model.RDS")
chrY_model <- read.csv("/path/to/sex_prediction_model/output/Yclass.predictions.csv", row.names = 1)
Xist_model <- read.csv("/path/to/sex_prediction_model/output/Xist.class.predictions.csv", row.names = 1)

### testing data
testing.VTA <- readRDS("/path/to/sex_prediction_model/data/Rn7_VTA_testing.RDS") #Seurat object
testing <- read.csv("/path/to/sex_prediction_model/data/test_count_data.csv", row.names = 1) #Read count table
#### make character ID bin
testing$Identity_bin_char <- recode(testing$Identity_bin, `0` = "Male", `1` = "Female")
testing$Identity_bin_char <- factor(testing$Identity_bin_char, levels = c("Male", "Female"))
#### make ID bin a factor
testing$Identity_bin <- as.factor(testing$Identity_bin)
#### add neuronal/nonneuronal celltype column
testing.VTA@meta.data <- testing.VTA@meta.data %>% 
  mutate(Neuronal = case_when(CellType %in% c("Glut-Neuron-1", "Glut-Neuron-2", "Glut-Neuron-3", "GABA-Neuron-1", "GABA-Neuron-2", "GABA-Neuron-3", "DA-Neuron") ~ "Neuronal",
                                  CellType %in% c("Olig-1", "Olig-2", "Olig-3", "Astrocyte", "Polydendrocyte", "Microglia", "OPC-Olig-1", "Mural", "Endothelial") ~ "Non-Neuronal"))

# Predictions ####
## logistic regression
lr_prediction <- predict(lr_model,testing[,1:285], type = "response")
## support vector machine
svm_prediction <- predict(svm_model,testing[,1:285], type = "prob")
## random forrest
rf_prediction <- predict(rf_model,testing[,1:285], type = "prob")
## multi-layer perceptron
mlp_prediction <- predict(mlp_model,testing[,1:285], type = "prob")
## chrY expression based prediction
chrY_prediction <- data.frame(cell = names(Idents(testing.VTA)),
                              sex = testing.VTA$Sex,
                              counts = y.count(testing.VTA),
                              prediction = y.count.classify(testing.VTA,0))
chrY_prediction$prediction <- factor(chrY_prediction$prediction, levels = c("Male", "Female"))
## Xist expression based prediction
Xist_prediction <- data.frame(cell = names(Idents(testing.VTA)),
                              sex = testing.VTA$Sex,
                              counts = Xist.count(testing.VTA),
                              prediction = Xist.count.classify(testing.VTA,0))
Xist_prediction$prediction <- factor(Xist_prediction$prediction, levels = c("Male", "Female"))

## make precitions using threshold
svm_prediction$class <- ifelse(svm_prediction$Female >0.5, "Female", "Male")
svm_prediction$class <- factor(svm_prediction$class, levels = c("Male", "Female"))

rf_prediction$class <- ifelse(rf_prediction$Female >0.5, "Female", "Male")
rf_prediction$class <- factor(rf_prediction$class, levels = c("Male", "Female"))

mlp_prediction$class <- ifelse(mlp_prediction$Female >0.5, "Female", "Male")
mlp_prediction$class <- factor(mlp_prediction$class, levels = c("Male", "Female"))

lr_prediction <- ifelse(lr_prediction >0.5, "Female", "Male")
lr_prediction <- factor(lr_prediction, levels = c("Male", "Female"))

# Aggregate results ####
acc.depth.df <- data.frame(cell = row.names(testing),
                           celltype = testing.VTA$Neuronal,
                           RNA = testing.VTA$nCount_RNA,
                           sex = testing.VTA$Sex,
                           chrY = chrY_prediction$prediction,
                           Xist = Xist_prediction$prediction,
                           lr = lr_prediction,
                           svm = svm_prediction$class,
                           rf = rf_prediction$class,
                           mlp = mlp_prediction$class)

## change predictions to correct?(T/F)
acc.depth.df <- acc.depth.df %>% 
  mutate(chrY = case_when(chrY == sex ~ 1,
                          chrY != sex ~ 0),
         Xist = case_when(Xist == sex ~ 1,
                          Xist != sex ~ 0),
         lr   = case_when(lr == sex ~ 1,
                          lr != sex ~ 0),
         svm  = case_when(svm == sex ~ 1,
                          svm != sex ~ 0),
         rf   = case_when(rf == sex ~ 1,
                          rf != sex ~ 0),
         mlp  = case_when(mlp == sex ~ 1,
                          mlp != sex ~ 0))

# Plot accuracy depth ####
## calculate logistic regressions for read depth and classification
fits <- list(chrY = glm(chrY ~ RNA, data = acc.depth.df, family = "binomial"),
             Xist = glm(Xist ~ RNA, data = acc.depth.df, family = "binomial"),
             lr = glm(lr ~ RNA, data = acc.depth.df, family = "binomial"),
             svm = glm(svm ~ RNA, data = acc.depth.df, family = "binomial"),
             rf = glm(rf ~ RNA, data = acc.depth.df, family = "binomial"),
             mlp = glm(mlp ~ RNA, data = acc.depth.df, family = "binomial"))

## collect results in data frame
model.intercept <- sapply(fits, function(x){coef(summary(x))[1,]}) %>% t() %>% data.frame()
model.intercept$model <- row.names(model.intercept)
colnames(model.intercept) <- c("Intercept_Estimate", "Intercept_StdErr", "Intercept_Z", "Interpect_pval", "model")

model.RNA <- sapply(fits, function(x){coef(summary(x))[2,]}) %>% t() %>% data.frame()
model.RNA$model <- row.names(model.RNA)
colnames(model.RNA) <- c("RNA_Estimate", "RNA_StdErr", "RNA_Z", "RNA_pval", "model")

results.df <- left_join(x = model.intercept,
                        y = model.RNA,
                        by = "model" ) %>% 
  select(model, Intercept_Estimate, Intercept_StdErr, Intercept_Z, Interpect_pval, RNA_Estimate, RNA_StdErr, RNA_Z, RNA_pval)

## calculate read depth at which probability of correct prediction is 95%
critical <- lapply(fits, function(model){
  ### Given values
  intercept <- model[["coefficients"]][["(Intercept)"]]
  coef_RNA <- model[["coefficients"]][["RNA"]]
  target_probability <- 0.95
  
  ### Define the equation to solve
  equation <- function(RNA) {
    p <- 1 / (1 + exp(-(intercept + coef_RNA * RNA)))
    return(p - target_probability)
  }
  
  ### Use uniroot to find the RNA value
  result <- uniroot(equation, c(0, 80000))  # Adjust the range as needed
  
  ### The result contains the RNA value
  return(round(result$root))
})

## make plots
theme_set(theme_bw())
theme_update(panel.grid.minor = element_blank(),
             panel.grid.major = element_blank(),
             panel.border = element_rect(colour = "black"),
             text = element_text(size = 10),
             axis.title = element_text(size = 8),
             axis.text = element_text(size = 7))

plots <- list(
  chrY_plot = ggplot(data = acc.depth.df,
                     mapping = aes(x = RNA, y = chrY)) +
    geom_point(shape = 1) +
    geom_smooth(method = "glm",
                method.args = list(family = "binomial"), se = F) +
    geom_hline(yintercept =  0.95, linetype = "dashed", color = "red") +
    geom_vline(xintercept =  critical$chrY, linetype = "dashed", color = "red") +
    scale_x_continuous(breaks = c(0, critical$chrY, 20000, 40000, 60000, 80000), guide = guide_axis(n.dodge = 2)) +
    scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 0.95, 1)) +
    labs(title = "Chr Y", x = "RNA count", y = "Classification"),
  
  Xist_plot = ggplot(data = acc.depth.df,
                     mapping = aes(x = RNA, y = Xist)) +
    geom_point(shape = 1) +
    geom_smooth(method = "glm",
                method.args = list(family = "binomial"), se = F) +
    geom_hline(yintercept =  0.95, linetype = "dashed", color = "red") +
    geom_vline(xintercept =  critical$Xist, linetype = "dashed", color = "red") +
    scale_x_continuous(breaks = c(0, critical$Xist, 20000, 40000, 60000, 80000),guide = guide_axis(n.dodge = 2)) +
    scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 0.95, 1)) +
    labs(title = "Xist", x = "RNA count", y = "Classification"),
  
  lr_plot = ggplot(data = acc.depth.df,
                   mapping = aes(x = RNA, y = lr)) +
    geom_point(shape = 1) +
    geom_smooth(method = "glm",
                method.args = list(family = "binomial"), se = F) +
    geom_hline(yintercept =  0.95, linetype = "dashed", color = "red") +
    geom_vline(xintercept =  critical$lr, linetype = "dashed", color = "red") +
    scale_x_continuous(breaks = c(0, critical$lr, 20000, 40000, 60000, 80000), guide = guide_axis(n.dodge = 2)) +
    scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 0.95, 1)) +
    labs(title = "LR", x = "RNA count", y = "Classification"),
  
  svm_plot = ggplot(data = acc.depth.df,
                    mapping = aes(x = RNA, y = svm)) +
    geom_point(shape = 1) +
    geom_smooth(method = "glm",
                method.args = list(family = "binomial"), se = F) +
    geom_hline(yintercept =  0.95, linetype = "dashed", color = "red") +
    geom_vline(xintercept =  critical$svm, linetype = "dashed", color = "red") +
    scale_x_continuous(breaks = c(0, critical$svm, 20000, 40000, 60000, 80000), guide = guide_axis(n.dodge = 2)) +
    scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 0.95, 1)) +
    labs(title = "SVM", x = "RNA count", y = "Classification"),
  rf_plot = ggplot(data = acc.depth.df,
                   mapping = aes(x = RNA, y = rf)) +
    geom_point(shape = 1) +
    geom_smooth(method = "glm",
                method.args = list(family = "binomial"), se = F) +
    geom_hline(yintercept =  0.95, linetype = "dashed", color = "red") +
    geom_vline(xintercept =  critical$rf, linetype = "dashed", color = "red") +
    scale_x_continuous(breaks = c(0, critical$rf, 20000, 40000, 60000, 80000), guide = guide_axis(n.dodge = 2)) +
    scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 0.95, 1)) +
    labs(title = "RF", x = "RNA count", y = "Classification"),
  mlp_plot = ggplot(data = acc.depth.df,
                    mapping = aes(x = RNA, y = mlp)) +
    geom_point(shape = 1) +
    geom_smooth(method = "glm",
                method.args = list(family = "binomial"), se = F) +
    geom_hline(yintercept =  0.95, linetype = "dashed", color = "red") +
    geom_vline(xintercept =  critical$mlp, linetype = "dashed", color = "red") +
    scale_x_continuous(breaks = c(0, critical$mlp, 20000, 40000, 60000, 80000), guide = guide_axis(n.dodge = 2)) +
    scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 0.95, 1)) +
    labs(title = "MLP", x = "RNA count", y = "Classification"))

## Arrange plots into figure
figure <- plots$chrY_plot + plots$Xist_plot + plots$lr_plot + plots$svm_plot + plots$rf_plot + plots$mlp_plot +
  plot_layout(ncol = 2)

# log scaled ####
plots.log <- list(
  chrY_plot = ggplot(data = acc.depth.df,
                     mapping = aes(x = RNA, y = chrY)) +
    geom_point(shape = 1) +
    geom_smooth(method = "glm",
                formula = y ~ I(10^x),
                method.args = list(family = "binomial"),
                se = F) + 
    geom_hline(yintercept =  0.95, linetype = "dashed", color = "red") +
    geom_vline(xintercept =  critical$chrY, linetype = "dashed", color = "red") +
    scale_x_continuous(breaks = c(0, 500, 1000,  5000, critical$chrY, 10000, 50000),
                       trans = "log10", guide = guide_axis(n.dodge = 2)) +
    scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 0.95, 1)) +
    labs(title = "Chr Y", x = "RNA count", y = "Classification"),
  
  Xist_plot = ggplot(data = acc.depth.df,
                     mapping = aes(x = RNA, y = Xist)) +
    geom_point(shape = 1) +
    geom_smooth(method = "glm",
                formula = y ~ I(10^x),
                method.args = list(family = "binomial"),
                se = F) +
    geom_hline(yintercept =  0.95, linetype = "dashed", color = "red") +
    geom_vline(xintercept =  critical$Xist, linetype = "dashed", color = "red") +
    scale_x_continuous(breaks = c(0, 500, 1000,  5000, 10000, critical$Xist, 50000),
                       trans = "log10", guide = guide_axis(n.dodge = 2)) +
    scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 0.95, 1)) +
    labs(title = "Xist", x = "RNA count", y = "Classification"),
  lr_plot = ggplot(data = acc.depth.df,
                   mapping = aes(x = RNA, y = lr)) +
    geom_point(shape = 1) +
    geom_smooth(method = "glm",
                formula = y ~ I(10^x),
                method.args = list(family = "binomial"),
                se = F) +
    geom_hline(yintercept =  0.95, linetype = "dashed", color = "red") +
    geom_vline(xintercept =  critical$lr, linetype = "dashed", color = "red") +
    scale_x_continuous(breaks = c(0, 500, 1000,  5000, critical$lr, 10000, 50000),
                       trans = "log10", guide = guide_axis(n.dodge = 2)) +
    scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 0.95, 1)) +
    labs(title = "LR", x = "RNA count", y = "Classification"),
  svm_plot = ggplot(data = acc.depth.df,
                    mapping = aes(x = RNA, y = svm)) +
    geom_point(shape = 1) +
    geom_smooth(method = "glm",
                formula = y ~ I(10^x),
                method.args = list(family = "binomial"),
                se = F) +
    geom_hline(yintercept =  0.95, linetype = "dashed", color = "red") +
    geom_vline(xintercept =  critical$svm, linetype = "dashed", color = "red") +
    scale_x_continuous(breaks = c(0, 500, 1000,  5000, critical$svm, 10000, 50000),
                       trans = "log10", guide = guide_axis(n.dodge = 2)) +
    scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 0.95, 1)) +
    labs(title = "SVM", x = "RNA count", y = "Classification"),
  rf_plot = ggplot(data = acc.depth.df,
                   mapping = aes(x = RNA, y = rf)) +
    geom_point(shape = 1) +
    geom_smooth(method = "glm",
                formula = y ~ I(10^x),
                method.args = list(family = "binomial"),
                se = F) +
    geom_hline(yintercept =  0.95, linetype = "dashed", color = "red") +
    geom_vline(xintercept =  critical$rf, linetype = "dashed", color = "red") +
    scale_x_continuous(breaks = c(0, 500, 1000,  5000, critical$rf, 10000, 50000),
                       trans = "log10", guide = guide_axis(n.dodge = 2)) +
    scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 0.95, 1)) +
    labs(title = "RF", x = "RNA count", y = "Classification"),
  mlp_plot = ggplot(data = acc.depth.df,
                    mapping = aes(x = RNA, y = mlp)) +
    geom_point(shape = 1) +
    geom_smooth(method = "glm",
                formula = y ~ I(10^x),
                method.args = list(family = "binomial"),
                se = F) +
    geom_hline(yintercept =  0.95, linetype = "dashed", color = "red") +
    geom_vline(xintercept =  critical$mlp, linetype = "dashed", color = "red") +
    scale_x_continuous(breaks = c(0, 500, 1000,  5000, critical$mlp, 10000, 50000),
                       trans = "log10", guide = guide_axis(n.dodge = 2)) +
    scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 0.95, 1)) +
    labs(title = "MLP", x = "RNA count", y = "Classification"))

## Arrange plots into figure
figure.log <- plots.log$chrY_plot + plots.log$Xist_plot + plots.log$lr_plot + plots.log$svm_plot + plots.log$rf_plot + plots.log$mlp_plot +
  plot_layout(ncol = 2)


#Save outputs
# logistic regression estimates
write.csv(results.df,
          file = "/path/to/sex_prediction_model/output/depth_vs_acc/depth_vs_acc.celltype.coeff.csv",
          row.names = F)

## non-scaled
pdf(file = "/path/to/sex_prediction_model/output/depth_vs_acc/depth_vs_acc.VTA.allmodels.pdf",
    width = 8,
    height = 8)
figure
dev.off()

## log10 scaled
pdf(file = "/path/to/sex_prediction_model/output/depth_vs_acc/depth_vs_acc.log.VTA.allmodels.pdf",
    width = 8,
    height = 8)
figure.log
dev.off()

# sessionInfo ####
sessionInfo()