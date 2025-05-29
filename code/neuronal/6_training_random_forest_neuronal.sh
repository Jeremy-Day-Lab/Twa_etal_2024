#!/bin/bash
#
#SBATCH --job-name=6_training_random_forest_neuronal
#SBATCH --output=/scratch/gtwa/Day/sex_prediction_model/Rn7_Xist_add/logs/neuronal/6_training_random_forest_neuronal.out
#SBATCH --error=/scratch/gtwa/Day/sex_prediction_model/Rn7_Xist_add/logs/neuronal/6_training_random_forest_neuronal.err
#SBATCH --ntasks=1
#SBATCH --partition=medium
#SBATCH --time=50:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=16g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=gtwa@uab.edu

# Load Modules
module load R/4.2.0-foss-2022b 

# Run script
Rscript /scratch/gtwa/Day/sex_prediction_model/Rn7_Xist_add/code/neuronal/6_training_random_forest_neuronal.R