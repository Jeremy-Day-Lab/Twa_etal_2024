# Accurate Sample Seconvolution of Pooled snRNA-seq Using Sex-Dependent Gene Expression Patterns
This repository contains the code for training, testing, and analysis of cell-sex prediction models reported in Twa et al., 2024. Any questions should be directed to jjday@uab.edu and gtwa@uab.edu.

## Citation
Guy M Twa, Robert A Phillips III, Nathaniel J Robinson, Jeremy J Day. bioRxiv 2024.11.29.626066; doi: https://doi.org/10.1101/2024.11.29.626066

## Raw Data
Sequencing data that support the findings of this study are available in Gene Expression Omnibus. Accession numbers of specific datasets are outlined below.
Ventral tegmental area snRNA-seq VTA: GSE168156
Nucleus accumbens snRNA-seq: GSE137763, GSE222418

## Models
Pre-trained models are available within the `models` directory. The models are trained on the rat VTA dataset as described in the manuscript. The models are saved in RDS format, they can be loaded in R using the `readRDS` function. The models can be used for cell sex prediction with a cell by gene expression matrix using the `predict()` function, as showin in code/9_evaluation_overall.R and other evaluation scripts.

If you download the models for your own use, please cite the original publication as described above. Thanks!