#! /bin/bash
dat=(BLCA COADREAD GBMLGG HNSC KIPAN LUAD THCA UCEC CESC KIRC KIRP LAML SARC ACC CHOL DLBC GBM LUSC READ THYM UCS UVM)
DATA="/icgc/dkfzlsdf/analysis/B080/steiger"

for data in ${dat[*]}; do
  echo "Plotting for $data"
  Rscript-3.1.2 --vanilla $DATA/scripts/Plot.R $data

done
