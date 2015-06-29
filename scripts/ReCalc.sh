#HighMem
dat=(BLCA COADREAD GBMLGG HNSC KIPAN LUAD THCA UCEC)
for data in ${dat[*]}; do
  qsub -v data="$data" -l\walltime=3:00:00,mem=45g CompLoad.pbs

done

#MediumMem
dat=(CESC KIRC KIRP LAML SARC)
for data in ${dat[*]}; do
  qsub -v data="$data" -l\walltime=3:00:00,mem=25g CompLoad.pbs

done

#LowMem
dat=(ACC CHOL DLBC GBM LUSC READ THYM UCS UVM)
for data in ${dat[*]}; do
  qsub -v data="$data" -l\walltime=3:00:00,mem=15g CompLoad.pbs

done
