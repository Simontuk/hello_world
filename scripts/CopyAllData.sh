#! /bin/bash
ssh steiger@b080-pc58 ls /data/TCGA > ALL

while read d; do

echo "Collecting $d TCGA Data \n"

DATA="/icgc/dkfzlsdf/analysis/B080/steiger"
cd "$DATA/data"

if [ -d $d ]; then
  continue
fi
mkdir $d
METH="gdac.broadinstitute.org_$d.Merge_methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.Level_3.2015020400.0.0.tar.gz"
RNASEQ="gdac.broadinstitute.org_$d.mRNAseq_Preprocess.Level_3.2015020400.0.0.tar.gz"

scp steiger@b080-pc58:"/data/TCGA/$d/stddata__2015_02_04/$d/20150204/$METH /data/TCGA/$d/stddata__2015_02_04/$d/20150204/$RNASEQ" $d/

cd $d
METHOUT="${METH%.tar.gz}/$d.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt"
RNASEQOUT="${RNASEQ%.tar.gz}/$d.uncv2.mRNAseq_RSEM_normalized_log2_PARADIGM.txt"
tar -zkxvf $METH $METHOUT
tar -zkxvf $RNASEQ -v --wildcards $RNASEQOUT

awk -F'\t' '{for(i=2;i<= NF; i+=4) if (i==2) printf $1"\t"; else printf ("%s%c", $i, i+4 <=NF ? "\t" :"\n");}' $METHOUT > methylation450_NArm.txt


done < ALL
