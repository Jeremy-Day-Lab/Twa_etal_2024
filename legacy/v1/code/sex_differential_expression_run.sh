#!/bin/bash
#
#SBATCH --job-name=sex_diff
#SBATCH --output=/path/to/sex_prediction_model/logs/sex_differential_expression_out.txt
#SBATCH --error=/path/to/sex_prediction_model/logs/sex_differential_expression_error.txt
#SBATCH --ntasks=1
#SBATCH --partition=short
#SBATCH --time=6:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=16g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=gtwa@uab.edu

# Load Modules
module load R/4.2.0-foss-2021a-bare

# Run script
Rscript /path/to/sex_prediction_model/sex-prediction-model/scripts/sex_differential_expression.R