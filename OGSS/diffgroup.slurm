#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=60
#SBATCH --time=00:01:00
#SBATCH --partition=msismall
#SBATCH --mem=100g
#SBATCH --account=tianxili
#SBATCH --output=out/datagen.out
#SBATCH -M agate


module load R
Rscript diffgroup.R