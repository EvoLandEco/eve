#!/bin/bash
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --job-name=edd_load
#SBATCH --output=logs/edd_load-%j.log
#SBATCH --mem=300GB
#SBATCH --partition=himem

name=${1}

ml R
Rscript ~/eve/Script/run_edd_load.R ${name}