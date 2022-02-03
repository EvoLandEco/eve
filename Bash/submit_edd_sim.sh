#!/bin/bash
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=edd_sim
#SBATCH --output=logs/edd_sim-%j.log
#SBATCH --mem=1GB
#SBATCH --partition=regular

name=${1}
param_set=${2}
la=${3}
mu=${4}
beta_n=${5}
beta_phi=${6}
gamma_n=${7}
gamma_phi=${8}
age=${9}
model=${10}
metric=${11}
offset=${12}

ml R
Rscript ~/eve/Script/run_edd_sim.R ${name} \
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
