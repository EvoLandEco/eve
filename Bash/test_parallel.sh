#!/bin/bash
#SBATCH --time=0:29:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --job-name=test_parallel
#SBATCH --output=logs/test_parallel-%j.log
#SBATCH --mem=10GB
#SBATCH --partition=short

ml R

Rscript ~/eve/Script/test_parallel.R
