#! /bin/bash
echo "Collecting $1 TCGA Data \n"

DATA="/icgc/dkfzlsdf/analysis/B080/steiger"
cd "$DATA/data"
mkdir $1
METH="gdac.broadinstitute.org_$1.Merge_methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.Level_3.2015020400.0.0.tar.gz"
RNASEQ="gdac.broadinstitute.org_$1.mRNAseq_Preprocess.Level_3.2015020400.0.0.tar.gz"

scp steiger@b080-pc58:"/data/TCGA/$1/stddata__2015_02_04/$1/20150204/$METH /data/TCGA/$1/stddata__2015_02_04/$1/20150204/$RNASEQ" $1/

cd $1
tar -zkxvf $METH "${METH%.tar.gz}/$1.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt" 
tar -zkxvf $RNASEQ -v --wildcards "${RNASEQ%.tar.gz}/$1.uncv2.mRNAseq_RSEM_normalized_log2_PARADIGM.txt"

perl "$DATA/scripts/ExtractMethylationDataNArm.pl" "$1"
#find . -type f -exec chmod a-w {} +


Rscript --vanilla $DATA/scripts/01_loadData.ALL.R $1


