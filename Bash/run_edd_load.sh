#!/bin/bash
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --job-name=edd_load
#SBATCH --output=logs/edd_load-%j.log
#SBATCH --mem=32GB
#SBATCH --partition=regular

name=${1}

ml R
Rscript ~/eve/Script/run_edd_load.R ${name}