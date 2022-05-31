#!/bin/bash
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --job-name=edd_analysis
#SBATCH --output=logs/edd_analysis-%j.log
#SBATCH --mem=80GB
#SBATCH --partition=regular



ml R

name=${1}

Rscript -e "devtools::install_github('EvoLandEco/eve')"

Rscript ~/eve/Script/run_edd_sim.R ${name}
