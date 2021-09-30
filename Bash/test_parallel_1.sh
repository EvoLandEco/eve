#!/bin/bash
#SBATCH --time=8:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --job-name=test_parallel
#SBATCH --output=logs/test_parallel-%j.log
#SBATCH --mem=10GB
#SBATCH --partition=regular

ml R

test_name=$1

Rscript ~/eve/Script/test_parallel.R ${test_name}
