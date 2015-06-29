#! /bin/bash

#while read data; do
# #echo $data
#
#done < /icgc/dkfzlsdf/analysis/B080/steiger/data/Cor.tmp
#! /bin/bash

#HighMem
dat=(BLCA COADREAD GBMLGG HNSC KIPAN LUAD THCA UCEC)
for data in ${dat[*]}; do
  qsub -v data="$data" -q highmem -l\walltime=200:00:00,mem=400g,nodes=1:ppn=12 CorCalc.pbs

done

#MediumMem
dat=(CESC KIRC KIRP LAML SARC)
for data in ${dat[*]}; do
  qsub -v data="$data" -l\walltime=200:00:00,mem=150g,nodes=1:ppn=12 CorCalc.pbs

done

#LowMem
dat=(ACC CHOL DLBC GBM LUSC READ THYM UCS UVM)
for data in ${dat[*]}; do
  qsub -v data="$data" -l\walltime=200:00:00,mem=80g,nodes=1:ppn=12 CorCalc.pbs

done
