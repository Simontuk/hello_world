#! /bin/bash
dat=(KIPAN OV READ SARC THCA THYM UCEC UVM)

for data in ${dat[*]}; do
  qsub -v data="$data" CompStat.pbs
  #qsub -v data="$data" -W depend=afterok:$FIRST CorCalc.pbs
done
