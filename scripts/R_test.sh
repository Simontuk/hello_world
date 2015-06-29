#!/bin/bash
#PBS -o R_out.log
#PBS -j oe

#PBS -d /icgc/dkfzlsdf/analysis/B080/steiger
#PBS -l walltime=1:00:00,mem=20g
#PBS -N R_test

echo We are in the $PWD directory
cd $PBS_O_WORKDIR
cd scripts/
echo We are now in $PWD, running an R script

R --vanilla < 01_loadData.ACC.R


