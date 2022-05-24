#!/bin/bash
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=edd_sim
#SBATCH --output=logs/edd_sim-%j.log
#SBATCH --mem=3GB
#SBATCH --partition=regular

name=${1}
param_set=${2}
nrep=${3}

ml R
Rscript ~/eve/Script/run_edd_sim.R ${name} \
                                   ${param_set} \
                                   ${nrep}
