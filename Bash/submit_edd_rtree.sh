#!/bin/bash
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=edd_sim
#SBATCH --output=logs/edd_sim-%j.log
#SBATCH --mem=4GB
#SBATCH --partition=regular

name=${1}
param_set=${2}
size=${3}
n=${4}

ml R
Rscript ~/eve/Script/run_edd_rtree.R ${name} \
                                   ${param_set} \
                                   ${size} \
                                   ${n}
