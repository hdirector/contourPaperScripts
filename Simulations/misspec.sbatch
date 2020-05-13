#!/bin/bash
#SBATCH --job-name=misspec
#SBATCH --partition=short
#SBATCH --time 0-11:59 # Runtime in D-HH:MM
#SBATCH --mem-per-cpu=50000
#SBATCH --out aMS_%j.out # File to which STDOUT will be written
#SBATCH --error aMS_%j.err # File to which STDERR will be written
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
#SBATCH --mail-user=direch@uw.edu
#SBATCH --array=1-6

module load R
Rscript runSimA.R $SLURM_ARRAY_TASK_ID "misspec" 

