#!/bin/bash
#SBATCH --time=00:10:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=edd_start
#SBATCH --output=logs/edd_start-%j.log
#SBATCH --mem=1GB
#SBATCH --partition=short

ml R

Rscript -e "devtools::install_github('EvoLandEco/eve')"

name=${1}
size=${2}
n=${3}

for (( param_set = 1; param_set <= 405; param_set++ ))
do
sbatch ~/eve/Bash/submit_edd_rtree.sh ${name} \
                                    ${param_set} \
                                    ${size} \
                                    ${n}
done
