#!/bin/bash
#SBATCH --job-name=nObsnGen
#SBATCH --partition=medium
#SBATCH --time 3-00:00 # Runtime in D-HH:MM
#SBATCH --mem-per-cpu=80000
#SBATCH --out OG_%j.out # File to which STDOUT will be written
#SBATCH --error OG_%j.err # File to which STDERR will be written
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
#SBATCH --mail-user=direch@uw.edu
#SBATCH --array=1-9

module load R
Rscript runSim.R $SLURM_ARRAY_TASK_ID "nObsnGen" 
