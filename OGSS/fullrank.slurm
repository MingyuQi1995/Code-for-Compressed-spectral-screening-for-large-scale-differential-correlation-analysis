#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --time=5:00:00
#SBATCH --partition=msismall
#SBATCH --mem=200g
#SBATCH --account=tianxili
#SBATCH --output=out/full.out
#SBATCH -M mesabi


module load R
Rscript fullrank.R