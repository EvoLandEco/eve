#!/bin/bash
#SBATCH --time=00:10:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=edd_start
#SBATCH --output=logs/edd_start-%j.log
#SBATCH --mem=1GB
#SBATCH --partition=short

ml R

Rscript -e "devtools::install_github('rsetienne/DDD@tianjian_Rampal')"
Rscript -e "devtools::install_github('EvoLandEco/eve')"

name=${1}
nrep=${2}

for (( param_set = 1; param_set <= 450; param_set++ ))
do
sbatch ~/eve/Bash/submit_edd_sim.sh ${name} \
                                    ${param_set} \
                                    ${nrep}
done
