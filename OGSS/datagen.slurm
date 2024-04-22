#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=26
#SBATCH --time=10:00:00
#SBATCH --partition=msismall,agsmall
#SBATCH --mem=200g
#SBATCH --account=tianxili
#SBATCH --output=out/datagen.out
#SBATCH -M agate


module load R
Rscript datagen.R