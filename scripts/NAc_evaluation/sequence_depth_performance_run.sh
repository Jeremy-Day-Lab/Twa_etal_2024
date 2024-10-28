#!/bin/bash
#
#SBATCH --job-name=depth_performance
#SBATCH --output=/scratch/gtwa/Day/sex_prediction_model/logs/depth_performance_out.txt
#SBATCH --error=/scratch/gtwa/Day/sex_prediction_model/logs/depth_performance_err.txt
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
Rscript /scratch/gtwa/Day/sex_prediction_model/sex-prediction-model/scripts/NAc_evaluation/sequence_depth_performance.R