#! /bin/bash

while read data; do
  echo $data
  #qsub -v data="$data" CorCalc.pbs
done < /icgc/dkfzlsdf/analysis/B080/steiger/data/Cor.tmp
