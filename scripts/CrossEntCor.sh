#PBS -m abe
#PBS -l walltime=250:00:00,mem=200g
#! /bin/bash

dat1=UCS
dat2=ACC

DATA="/icgc/dkfzlsdf/analysis/B080/steiger"
#cd "$DATA/data"

echo "Calcualting Correlation"
Rscript-3.1.2 --vanilla $DATA/scripts/CorMX2Ent.R $dat1 $dat2
