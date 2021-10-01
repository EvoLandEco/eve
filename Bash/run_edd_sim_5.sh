#!/bin/bash
#SBATCH --time=9-23:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --job-name=test_parallel
#SBATCH --output=logs/test_parallel-%j.log
#SBATCH --mem=10GB
#SBATCH --partition=regular

ml R

test_name=$1

Rscript ~/eve/Script/run_edd_sim.R ${test_name}
