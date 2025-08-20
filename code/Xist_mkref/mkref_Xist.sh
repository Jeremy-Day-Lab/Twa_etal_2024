#!/bin/bash
#
#SBATCH --job-name=Xist_ref
#SBATCH --output=Xist_ref.out
#SBATCH --error=Xist_ref.err
#SBATCH --ntasks=1
#SBATCH --partition=express
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=64g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=gtwa@uab.edu

# move to working directory
cd /data/project/daylab/Genome-Reference/Genomes/mRatBN7.2/Ensembl/release_105

# Remove overlapping gene 
awk '!/gene_id "ENSRNOG00000037911"/' Rattus_norvegicus.mRatBN7.2.105.gtf > Rattus_norvegicus.mRatBN7.2.105.Xist.gtf

# Append new Xist annotation
cat xist_annotation.gtf >> Rattus_norvegicus.mRatBN7.2.105.Xist.gtf

# Sort by chromosome and start position
sort -k1,1 -k4,4n Rattus_norvegicus.mRatBN7.2.105.Xist.gtf -o Rattus_norvegicus.mRatBN7.2.105.Xist.gtf

#Filter GTF
/home/gtwa/tools/cellranger/cellranger-6.1.2/cellranger mkgtf \
Rattus_norvegicus.mRatBN7.2.105.Xist.gtf Rattus_norvegicus.mRatBN7.2.105.Xist.filtered.gtf \
--attribute=gene_biotype:protein_coding \
--attribute=gene_biotype:lncRNA \
--attribute=gene_biotype:antisense \
--attribute=gene_biotype:IG_LV_gene \
--attribute=gene_biotype:IG_V_gene \
--attribute=gene_biotype:IG_V_pseudogene \
--attribute=gene_biotype:IG_D_gene \
--attribute=gene_biotype:IG_J_gene \
--attribute=gene_biotype:IG_J_pseudogene \
--attribute=gene_biotype:IG_C_gene \
--attribute=gene_biotype:IG_C_pseudogene \
--attribute=gene_biotype:TR_V_gene \
--attribute=gene_biotype:TR_V_pseudogene \
--attribute=gene_biotype:TR_D_gene \
--attribute=gene_biotype:TR_J_gene \
--attribute=gene_biotype:TR_J_pseudogene \
--attribute=gene_biotype:TR_C_gene

#Run mkref
/home/gtwa/tools/cellranger/cellranger-6.1.2/cellranger mkref \
--genome=mRatBN7_Xist \
--fasta=Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa \
--genes=Rattus_norvegicus.mRatBN7.2.105.Xist.filtered.gtf \
--ref-version=1.0.0

