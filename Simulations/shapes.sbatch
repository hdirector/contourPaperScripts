#!/bin/bash
#SBATCH --job-name=shapes
#SBATCH --partition=medium
#SBATCH --time 3-00:00 # Runtime in D-HH:MM
#SBATCH --mem-per-cpu=60000
#SBATCH --out SH_%j.out # File to which STDOUT will be written
#SBATCH --error SH_%j.err # File to which STDERR will be written
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
#SBATCH --mail-user=direch@uw.edu
#SBATCH --array=1-3

module load R
Rscript runSim.R $SLURM_ARRAY_TASK_ID "shapes" 
