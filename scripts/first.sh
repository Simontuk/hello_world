cd ../data/KIRP/gdac.broadinstitute.org_KIRP.Merge_methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.Level_3.2015020400.0.0/
awk -F'\t' '{for(i=2;i<= NF; i+=4) if (i==2) printf $1"\t"; else printf ("%s%c", $i, i+4 <=NF ? "\t" :"\n");}' KIRP.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt > ../methylation450_NArm.txt