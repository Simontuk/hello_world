#!/usr/bin/perl
use strict;
use warnings;

my $data = "/icgc/dkfzlsdf/analysis/B080/steiger/data/"; 
my $cancer =shift @ARGV or die "Usage $0 FILENAME\n";

my $input = $data."$cancer/gdac.broadinstitute.org_$cancer.Merge_methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.Level_3.2015020400.0.0/$cancer.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt";
my $output = $data."$cancer/methylation450_NArm.txt";

open (OUT,'>'.$output) or die "Could not open $output for writing: $!\n";
close(OUT);

open(TEXT, $input) or die "File $input  not found\n";


my @row = <TEXT>;
#
# for (my $ i = 0; $i <= 15; $i++){
#   my @line = split("\t",$col[$i]);
#   print $line[0]."\t".$line[1]."\t\t".$line[5]."\n";
# }
my @col;
my @tab;
my @head;


for (my $j = 0; $j<@row; $j++){
  @col = split("\t",$row[$j]);
  push @tab, [@col];
}
close (TEXT);

open (OUT,'>'.$output) or die "Could not open $output for writing: $!\n";

for (my $i = 0; $i<@row;$i++){
	if ($tab[$i]->[0] ne "Composite Element REF" && $tab[$i]->[1] ne "NA"){
	print OUT $tab[$i]->[0];
  	for (my $v = 0; $v<@col;$v++){
    		if ($tab[1]->[$v] eq "Beta_value" && $tab[$i]->[0] ne "Composite Element REF" && $tab[$i]->[$v] ne "NA"){
      		print OUT "\t".$tab[$i]->[$v];
    		}
  	}	
  print OUT "\n";
  }


}

close OUT;


# %%%%%%%%%%%%%%%
# Jeden Proband nur einmalig in @Head
# %%%%%%%%%%%%%%%
#
# for (my $k = 0; $k<10;$k++){
#     if ($k<=1){
#       push @head, $tab[0]->[$k];
#     }
#     elsif($head[-1] eq $tab[0]->[$k]){
#     }
#     else{
#       push @head, $tab[0]->[$k];
#     }
#
# }
# %%%%%%%%%%%%%%%%

# print join "\t",@head;
# print "\n";
    # if ($i <= 4){
    #   print $line->[$v]."\t";
    # }
    # else{
      # print "\n";
      # if ($line eq $tab[$v]->[0]){
      #
      #   print $tab[$v]->[0]."\t";
      #   $i = 0;
      # }
      # else{
      #   print $tab[$v]->[0]."\t".$line."\t";
      #   $i = 1;
      # }


    # }
