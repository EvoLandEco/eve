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
la=${2}
mu=${3}
beta_n=${4}
beta_phi=${5}
gamma_n=${6}
gamma_phi=${7}
age=${8}
model=${9}
metric=${10}
offset=${11}

nrow=`wc -l ~/eve/Data/${name}.csv | cut -f1 -d' '`
nrow=$(( ${nrow} - 1 ))

for (( param_set = 1; param_set <= $nrow; param_set++ ))
do
sbatch ~/eve/Bash/submit_edd_sim.sh ${name} \
                                    ${param_set} \
                                    ${la} \
                                    ${mu} \
                                    ${beta_n} \
                                    ${beta_phi} \
                                    ${gamma_n} \
                                    ${gamma_phi} \
                                    ${age} \
                                    ${model} \
                                    ${metric} \
                                    ${offset}
done
