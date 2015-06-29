#HighMem
dat1=(BLCA COADREAD GBMLGG HNSC KIPAN LUAD THCA UCEC CESC KIRC KIRP LAML SARC ACC CHOL DLBC GBM LUSC READ THYM UCS UVM)
dat2=(BLCA COADREAD GBMLGG HNSC KIPAN LUAD THCA UCEC CESC KIRC KIRP LAML SARC ACC CHOL DLBC GBM LUSC READ THYM UCS UVM)
x=0
cedat="/icgc/dkfzlsdf/analysis/B080/steiger/analysis/Cross-Entity"

for entry in "$cedat"/*
do
  datlist=($datlist"$(basename $entry)")
  s=$(($s+1))
done
#echo ${datlist[*]}

for data1 in ${dat1[*]}; do
  for data2 in ${dat2[*]}; do
    if [[ "$data1" != "$data2" ]]; then
      #echo "$data1""_""$data2"
      if [[ ("${dats[*]}" == *"$data1""_""$data2"*) || ("${dats[*]}" == *"$data2""_""$data1"*) ]];
      then
      #echo "False Pair"
      h=$(($h+1))
        if [[ ("${datlist[*]}" == *"$data1""_""$data2"*) || ("${datlist[*]}" == *"$data2""_""$data1"*)]]; then
          #echo "False Pair + Calced"
          v=$(($v+1))
          #echo $v
        fi
        #echo $dats

      else
        #echo "Right Pair"
        x=$(($x+1))
        #echo $x

        if [[ ("${datlist[*]}" == *"$data1""_""$data2"*) || ("${datlist[*]}" == *"$data2""_""$data1"*) ]]; then
          #echo "Right Pair + Not Calced"
          i=$(($i+1))
          #echo $i
        else
          echo "$data1""_""$data2"
          qsub -v data1="$data1",data2="$data2" -l\walltime=15:00:00,mem=60g,nodes=1:ppn=12 CECor.pbs
          j=$(($j+1))
        fi

        #
      fi
      dats=($dats"$data1""_""$data2")

      #echo "____"
    fi

  done

done
#
# echo "Sum of Cross-Entities Calculated"
# echo $s
#
# echo "False Pairs"
# echo $h
#
# echo "False Pairs + Calced"
# echo $v
#
# echo "Right Pairs"
# echo $x
#
# echo "Right Pair + Calced"
# echo $i
#
# echo "Right Pair + Not Calced"
# echo $j
