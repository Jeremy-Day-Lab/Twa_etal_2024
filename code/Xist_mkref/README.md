# INFO
This directory contains the code used to create the custom Ensembl mRatBn7.2 (Rn7, v105) genome reference containing annotation data for *Xist*, from the NCBI RefSeq annotation file (GCF_015227675.2).

## Ensembl mRatBN7.2 reference files
The Ensembl mRatBN7.2 reference files were obtained from the Ensembl FTP site at these locations:
GTF: https://ftp.ensembl.org/pub/release-105/gtf/rattus_norvegicus/Rattus_norvegicus.mRatBN7.2.105.gtf.gz
FASTA: https://ftp.ensembl.org/pub/release-105/fasta/rattus_norvegicus/dna/Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa.gz

md5sum files for the uncompressed reference files are provided in this directory as `mRatBN7_FASTA.md5` and `mRatBN7_GTF.md5`.

## Obtaining Xist gene annotation data
The Xist gene annotation data was obtained from the NCBI RefSeq database at this location:
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/015/227/675/GCF_015227675.2_mRatBN7.2/

Corresponding annotation lines for Xist were extracted from the GTF file using `grep "Xist" GCF_015227675.2_mRatBN7.2_genomic.gtf`. These lines were then manually edited to match the Ensembl formatting of `Rattus_norvegicus.mRatBN7.2.105.gtf` and saved as `xist_annotation.gtf`.

## Creating the custom genome annotation file
Genome annotations overlapping the *Xist* loci were manually identified and removed to gaurd against any potential issues in quantifying feature counts due to multi-mapping. The `xist_annotation.gtf` was then appended to the Ensembl GTF file and sorted by chromosome and start position.

## Cell Ranger mkref
The custom genome reference was then created following 10x Genomics custom reference creation guidelines.

## Feature file
To facilitate mapping of feature names from this analysis for others using Ensembl annotations, this directory also contains `VTA_F1_features.tsv.gz`. This is the cellranger count `features.tsv.gz` out file from analysis of the F1 GEM well of the VTA dataset. This file contains the Ensembl gene IDs, and gene names for all features detected.