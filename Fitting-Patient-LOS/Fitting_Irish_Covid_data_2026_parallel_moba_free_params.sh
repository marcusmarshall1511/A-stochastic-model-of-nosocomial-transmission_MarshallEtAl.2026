#!/bin/bash

# -- Slurm Directives --
#SBATCH --job-name=Fitting_Irish_Covid_data_2026_parallel_moba       # Job name
#SBATCH --output=logs/out_%A_%a.txt # Standard output and error log, %j is the job ID
#SBATCH --error=logs/err_%A_%a.txt   # Separate error log
#SBATCH --partition=nodes              
#SBATCH --time=24:00:00                  # Total run time limit (HH:MM:SS)
#SBATCH --ntasks=1                       # Run a single task
#SBATCH --cpus-per-task=10               # Number of CPU cores per task
#SBATCH --mem=3G                         # Job memory request (e.g., 4 GB)


# Load the Conda module
module load miniforge

# Activate your R environment
conda activate myRenv

module purge
module load R


# Execute your R script
echo "Running the R script..."
Rscript Fitting_Irish_Covid_data_2026_parallel_moba_free_params.R 