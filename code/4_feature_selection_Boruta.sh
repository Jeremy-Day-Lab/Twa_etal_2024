#!/bin/bash
#
#SBATCH --job-name=4_feature_selection_Boruta
#SBATCH --output=/scratch/gtwa/Day/sex_prediction_model/Rn7_Xist_add/logs/4_feature_selection_Boruta.out
#SBATCH --error=/scratch/gtwa/Day/sex_prediction_model/Rn7_Xist_add/logs/4_feature_selection_Boruta.err
#SBATCH --ntasks=1
#SBATCH --partition=long
#SBATCH --time=150:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=16g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=gtwa@uab.edu

# Load Modules
module load R/4.2.0-foss-2022b 

# Run script
Rscript /scratch/gtwa/Day/sex_prediction_model/Rn7_Xist_add/code/4_feature_selection_Boruta.R