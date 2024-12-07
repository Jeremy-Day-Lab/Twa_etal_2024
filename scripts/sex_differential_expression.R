# Setup ####
## libraries
library(Seurat)
#library(rtracklayer) #required, but we will only be using the import function
library(tidyverse)

## seed
set.seed(1234)

## data
### training object
Rn7_VTA <- readRDS("/path/to/sex_prediction_model/data/Rn7_VTA_training.RDS")
### Rn7 gene names
Rn7_gtf <- rtracklayer::import("/path/to/sex_prediction_model/data/Rattus_norvegicus.mRatBN7.2.105.gtf")
Rn7_gtf <- as.data.frame(Rn7_gtf)

# Differential Expression ####
## data setup
### set identities for cells
Idents(Rn7_VTA) <- Rn7_VTA$CellType_Sex

### initialize results list
Sex_DEGs <- vector(mode = "list",length = 16)
names(Sex_DEGs) <- levels(Rn7_VTA$CellType)
for(i in names(Sex_DEGs)){
  Sex_DEGs[[i]] <- FindMarkers(object = Rn7_VTA,
                               ident.1 = paste0(i,"_Female"),
                               ident.2 = paste0(i,"_Male"))
}
### flatten list
Sex_DEGs <- do.call(what = rbind,Sex_DEGs)

## Create an ID column that is the rownames
Sex_DEGs$ID <- rownames(Sex_DEGs)

## Add gene names and cluster
Sex_DEGs <- separate(Sex_DEGs, col = ID, into = c("Cluster", "GeneName"), sep = "\\.", remove = FALSE, extra = "merge")

### subset for significant DEGs
Sex_DEGs.significant <- subset(Sex_DEGs,subset=(p_val_adj < 0.05))

## Build a table of the frequency by which the gene is differentially expressed by sex
Sex_DEGs_Freq <- as.data.frame(table(Sex_DEGs.significant$GeneName))

## Get chromosomes of all Sex_DEGs
for(i in 1:nrow(Sex_DEGs_Freq)){
  if(is.na(Rn7_gtf[which(Rn7_gtf$gene_name  %in% as.character(Sex_DEGs_Freq[i,"Var1"]))[1],"gene_name"]) == TRUE){
    Sex_DEGs_Freq[i,"Chr"] <- as.character(Rn7_gtf[which(Rn7_gtf$gene_id  %in% as.character(Sex_DEGs_Freq[i,"Var1"]))[1],
                                                   "seqnames"])
  }else{
    Sex_DEGs_Freq[i,"Chr"] <- Rn7_gtf[which(Rn7_gtf$gene_name  %in% as.character(Sex_DEGs_Freq[i,"Var1"]))[1],"seqnames"]
  }
}

## Get the mean fold change of the genes
Sex_DEGs_Freq$Mean_logFC <- NA
for(i in 1:nrow(Sex_DEGs_Freq)){
  if(Sex_DEGs_Freq[i,"Freq"]>1){
    Sex_DEGs_Freq[i,"Mean_logFC"] <- mean(subset(Sex_DEGs.significant,subset=(GeneName == as.character(Sex_DEGs_Freq[i,"Var1"])))$avg_log2FC)
  }else{
    Sex_DEGs_Freq[i,"Mean_logFC"] <- subset(Sex_DEGs.significant,subset=(GeneName == as.character(Sex_DEGs_Freq[i,"Var1"])))$avg_log2FC
  }
}

## Pull counts matrix
### First pull the count data
count_data_Full <- as.data.frame(t(as.matrix(GetAssayData(object = Rn7_VTA,slot = "data",assay = "RNA"))))
### check if rownames are in same order as the identities vector 
all(names(Idents(Rn7_VTA)) == rownames(count_data_Full)) #TRUE 
### Now just create a new column in the count_data_full for identity
count_data_Full$Identity <- as.factor(Rn7_VTA$Sex)
### Convert identitity to binary code. 
count_data_Full$Identity_bin <- ifelse(count_data_Full$Identity  == "Female",
                                       1, #Females are 1
                                       0) #Males are 0



## Create a subset of count data full containing all of the Sex_DEGs
count_data_subset <- count_data_Full[,c(as.character(Sex_DEGs_Freq$Var1),"Identity_bin")]





# Save outputs ####
write.csv(Sex_DEGs, file = "/path/to/sex_prediction_model/data/DEGs_sex_cluster_VTA_full.csv", row.names = F)
write.csv(Sex_DEGs.significant, file = "/path/to/sex_prediction_model/data/DEGs_sex_cluster_VTA_sig.csv", row.names = F)


write.table(x         = Sex_DEGs_Freq,
            file      = "/path/to/sex_prediction_model/data/DEGsbySexbyCluster_VTA.txt",
            sep       = "\t",
            col.names = TRUE,
            row.names = FALSE,
            quote     = FALSE)

saveRDS(object = count_data_subset,file = "/path/to/sex_prediction_model/data/count_data_subset_100423.RDS")

# Session info ####
sessionInfo()