#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=20:00:00
#SBATCH --partition=ag2tb,msismall,agsmall,preempt
#SBATCH --mem=490g
#SBATCH --account=tianxili
#SBATCH --output=out/eval.out
#SBATCH -M agate


module load R
Rscript eval.R