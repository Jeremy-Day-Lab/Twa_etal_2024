#!/bin/bash
#
#SBATCH --job-name=15_depth_accuracy
#SBATCH --output=/scratch/gtwa/Day/sex_prediction_model/Rn7_Xist_add/logs/15_depth_accuracy.out
#SBATCH --error=/scratch/gtwa/Day/sex_prediction_model/Rn7_Xist_add/logs/15_depth_accuracy.err
#SBATCH --ntasks=1
#SBATCH --partition=express
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=16g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=gtwa@uab.edu

# Load Modules
module load R/4.2.0-foss-2022b 

# Run script
Rscript /scratch/gtwa/Day/sex_prediction_model/Rn7_Xist_add/code/15_depth_accuracy.R