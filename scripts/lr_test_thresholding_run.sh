#!/bin/bash
#
#SBATCH --job-name=lr_test_prediction_threshold
#SBATCH --output=/path/to/sex_prediction_model/logs/lr_test_prediction_threshold_out.txt
#SBATCH --error=/path/to/sex_prediction_model/logs/lr_test_prediction_threshold_error.txt
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
Rscript /path/to/sex_prediction_model/sex-prediction-model/scripts/lr_test_thresholding.R