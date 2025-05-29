#!/bin/bash
#
#SBATCH --job-name=feature_selection_plot
#SBATCH --output=/path/to/sex_prediction_model/logs/feature_selection_plot_out.txt
#SBATCH --error=/path/to/sex_prediction_model/logs/feature_selection_plot_error.txt
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
Rscript /path/to/sex_prediction_model/sex-prediction-model/scripts/feature_selection_plot.R