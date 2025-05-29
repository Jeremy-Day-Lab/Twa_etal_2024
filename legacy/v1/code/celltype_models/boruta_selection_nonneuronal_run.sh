#!/bin/bash
#
#SBATCH --job-name=boruta_selection_neuronal
#SBATCH --output=/scratch/gtwa/Day/sex_prediction_model/logs/boruta_selection_nonneuronal_out.txt
#SBATCH --error=/scratch/gtwa/Day/sex_prediction_model/logs/boruta_selection_nonneuronal_error.txt
#SBATCH --ntasks=1
#SBATCH --partition=medium
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=16g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=gtwa@uab.edu

# Load Modules
module load R/4.2.0-foss-2021a-bare

# Run script
Rscript /scratch/gtwa/Day/sex_prediction_model/sex-prediction-model/scripts/celltype_models/boruta_selection_nonneuronal.R