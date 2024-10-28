# Setup ####
## libraries
library(dplyr)
library(ggplot2)
library(tidyr)
library(forcats)

## data
# all.minfo.df <- read.csv("/path/to/sex_prediction_model/data/mutualinfo.df.csv", row.names = 1) %>%
#   mutate(cell.type = recode(cell.type, neuronal = "Neuronal", nonneuronal = "Non-neuronal")) %>% 
#   mutate(cell.type = factor(cell.type, levels = c("Neuronal", "Non-neuronal")))
training.minfo.df <- read.csv("/path/to/sex_prediction_model/data/mutualinfo.training.df.csv", row.names = 1) %>%
  mutate(cell.type = recode(cell.type, neuronal = "Neuronal", nonneuronal = "Non-neuronal")) %>% 
  mutate(cell.type = factor(cell.type, levels = c("Neuronal", "Non-neuronal")))
testing.minfo.df <- read.csv("/path/to/sex_prediction_model/data/mutualinfo.testing.df.csv", row.names = 1) %>%
  mutate(cell.type = recode(cell.type, neuronal = "Neuronal", nonneuronal = "Non-neuronal")) %>%
  mutate(cell.type = factor(cell.type, levels = c("Neuronal", "Non-neuronal")))

# Plot mutual info distributions ####
## histogram
training.hist <- ggplot(data = training.minfo.df, aes(x = mutual.info, fill = cell.type)) +
  geom_histogram(color = "black", bins = 30) +
  facet_wrap(~cell.type) +
  scale_x_continuous(trans = "log10") +
  labs(x = "Mutual Information", y = "Count") +
  theme_bw() +
  theme(legend.position = "none", text = element_text(size = 15))
testing.hist <- ggplot(data = testing.minfo.df, aes(x = mutual.info, fill = cell.type)) +
  geom_histogram(color = "black", bins = 30) +
  facet_wrap(~cell.type) +
  scale_x_continuous(trans = "log10") +
  labs(x = "Mutual Information", y = "Count") +
  theme_bw() +
  theme(legend.position = "none", text = element_text(size = 15))

### density
training.dens <- ggplot(data = training.minfo.df, aes(x = mutual.info, fill = cell.type)) +
  geom_density(alpha = 0.7) +
  scale_x_continuous(trans = "log10") +
  labs(x = "Mutual Information", y = "Density", fill = "Cell Type") +
  theme_bw() +
  theme(text = element_text(size = 15), legend.position = c(0.1,0.9))
testing.dens <- ggplot(data = testing.minfo.df, aes(x = mutual.info, fill = cell.type)) +
  geom_density(alpha = 0.7) +
  scale_x_continuous(trans = "log10") +
  labs(x = "Mutual Information", y = "Density", fill = "Cell Type") +
  theme_bw() +
  theme(text = element_text(size = 15), legend.position = c(0.1,0.9))

### violin
training.vln <- ggplot(data = training.minfo.df, aes(x = cell.type, y = mutual.info, fill = cell.type)) +
  geom_violin() +
  scale_y_continuous(trans = "log10") +
  labs(x = "Cell Type", y = "Mutual Information") +
  theme_bw() +
  theme(legend.position = "none", text = element_text(size = 15),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())
testing.vln <- ggplot(data = testing.minfo.df, aes(x = cell.type, y = mutual.info, fill = cell.type)) +
  geom_violin() +
  scale_y_continuous(trans = "log10") +
  labs(x = "Cell Type", y = "Mutual Information") +
  theme_bw() +
  theme(legend.position = "none", text = element_text(size = 15),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())

# statistical test ####
## KS test
### are these distributions significantly different than one another?
ks.training <- ks.test(mutual.info ~ cell.type, data = training.minfo.df)
ks.testing <- ks.test(mutual.info ~ cell.type, data = testing.minfo.df)

## Wilcox
### are the mutual information scores for genes higher in neuronal 
wilcox.training <- wilcox.test(mutual.info ~ cell.type, data = training.minfo.df, alternative = "g")
wilcox.testing <- wilcox.test(mutual.info ~ cell.type, data = testing.minfo.df, alternative = "g")


# paired compare
training.minfo.df.wide <- training.minfo.df %>% 
  pivot_wider(names_from = "cell.type",
              values_from = "mutual.info") %>% mutate(difference = Neuronal - `Non-neuronal`)

training.diff <- ggplot(training.minfo.df.wide, aes(x = fct_reorder(Gene, difference), y = difference)) +
  geom_point(shape = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Genes", y = "Mutual information (Neuronal - Non-neuronal)") +
  theme_bw() +
  theme(axis.text.x=element_blank(),
        panel.grid = element_blank())
  
testing.minfo.df.wide <- testing.minfo.df %>% 
  pivot_wider(names_from = "cell.type",
              values_from = "mutual.info") %>% mutate(difference = Neuronal - `Non-neuronal`)

testing.diff <- ggplot(testing.minfo.df.wide, aes(x = fct_reorder(Gene, difference), y = difference)) +
  geom_point(shape = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Genes", y = "Mutual information (Neuronal - Non-neuronal)") +
  theme_bw() +
  theme(axis.text.x=element_blank(),
        panel.grid = element_blank())

# write outputs ####
## distribution plots
### training
ggsave(plot = training.hist,
       file = "/path/to/sex_prediction_model/output/mutual_information/training.hist.pdf",
       device = "pdf",
       height = 10,
       width = 10)
ggsave(plot = training.dens,
       file = "/path/to/sex_prediction_model/output/mutual_information/training.dens.pdf",
       device = "pdf",
       height = 10,
       width = 10)
ggsave(plot = training.vln,
       file = "/path/to/sex_prediction_model/output/mutual_information/training.vln.pdf",
       device = "pdf",
       height = 10,
       width = 10)
ggsave(plot = training.diff,
       file = "/path/to/sex_prediction_model/output/mutual_information/training.diff.rank.pdf",
       device = "pdf",
       height = 10,
       width = 10)
write.csv(training.minfo.df.wide, file = "/path/to/sex_prediction_model/output/mutual_information/training.minfo.wide.csv")


### testing
ggsave(plot = testing.hist,
       file = "/path/to/sex_prediction_model/output/mutual_information/testing.hist.pdf",
       device = "pdf",
       height = 10,
       width = 10)
ggsave(plot = testing.dens,
       file = "/path/to/sex_prediction_model/output/mutual_information/testing.dens.pdf",
       device = "pdf",
       height = 10,
       width = 10)
ggsave(plot = testing.vln,
       file = "/path/to/sex_prediction_model/output/mutual_information/testing.vln.pdf",
       device = "pdf",
       height = 10,
       width = 10)
ggsave(plot = testing.diff,
       file = "/path/to/sex_prediction_model/output/mutual_information/testing.diff.rank.pdf",
       device = "pdf",
       height = 10,
       width = 10)
write.csv(testing.minfo.df.wide, file = "/path/to/sex_prediction_model/output/mutual_information/testing.minfo.wide.csv")

## print statistical tests
### training
cat("KS Test Training data: \n")
ks.training

cat("Wilcox Test Training data: \n")
wilcox.training

### testing
cat("KS Test Testing ata: \n")
ks.testing

cat("Wilcox Test Testing data: \n")
wilcox.testing


# sessionInfo ####
sessionInfo()