#!/bin/bash
#
#SBATCH --job-name=9_evaluation_overall_nonneuronal
#SBATCH --output=/scratch/gtwa/Day/sex_prediction_model/Rn7_Xist_add/logs/non_neuronal/9_evaluation_overall_nonneuronal.out
#SBATCH --error=/scratch/gtwa/Day/sex_prediction_model/Rn7_Xist_add/logs/non_neuronal/9_evaluation_overall_nonneuronal.err
#SBATCH --ntasks=1
#SBATCH --partition=express
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=16g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=user@mail.edu

# Load Modules
module load R/4.2.0-foss-2022b 

# Run script
Rscript /scratch/gtwa/Day/sex_prediction_model/Rn7_Xist_add/code/non_neuronal/9_evaluation_overall_nonneuronal.R