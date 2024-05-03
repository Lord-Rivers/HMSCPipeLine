#!/bin/bash
#SBATCH --job-name=KFold_Post_Processing
#SBATCH  -M ukko
#SBATCH --partition=short 
#SBATCH -c 2
#SBATCH --mem=32G
#SBATCH --time=00:10:00
#SBATCH --output=R-%x.%j.out

module purge
module load R

srun Rscript --vanilla --verbose S4_HPC_Post_Processing.R