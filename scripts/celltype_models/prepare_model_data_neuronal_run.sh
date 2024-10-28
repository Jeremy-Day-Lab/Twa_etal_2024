#!/bin/bash
#
#SBATCH --job-name=prepare_model_data_neuronal
#SBATCH --output=/scratch/gtwa/Day/sex_prediction_model/logs/prepare_model_data_neuronal_out.txt
#SBATCH --error=/scratch/gtwa/Day/sex_prediction_model/logs/prepare_model_data_neuronal_error.txt
#SBATCH --ntasks=1
#SBATCH --partition=express
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=16g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=gtwa@uab.edu

# Load Modules
module load R/4.2.0-foss-2021a-bare

# Run script
Rscript /scratch/gtwa/Day/sex_prediction_model/sex-prediction-model/scripts/celltype_models/prepare_model_data_neuronal.R