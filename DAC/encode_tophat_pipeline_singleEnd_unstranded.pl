#!/usr/bin/perl
#
# encode_tophat_pipeline.pl
# xwei 09/15/2014

## NOTE:: using tophat version -- tophat/2.0.8
##                                bowtie2/2.1.0
##                                samtools/0.1.17

## Modification 09/15/2014
# 1. for single-end and unstranded data
#       - remove "--no-discordant" and "--no-mixed" from tophat command, because they are for paired reads
#       - only have one input fastq file, $FASTQ1, for tophat
#       - change the $library_type to "fr-unstranded", for tophat
#       - no tophat_bam_xsA_tag_fix.pl needed because it works for pair-end stranded reads only
#       - using "unstranded" instead of "plus" or "minus" for parameters of encode_bam_to_bigwig.pl 

## Modification 04/25/2014 
# 1. using the agreed tophat parameters
#       --min-intron-length 20
#       --max-intron-length 1000000
#       -a 8
#       -m 0
#       remove "-n 2"
#       --read-edit-dist 4
#       --read-mismatches 4
# 2. edit the header of different bam files  
# 3. add bam to bigwig steps together

## Modification 09/17/2013
# 1. generate the index for accepted_hits.bam file
# 2. increase the running core number to 16

## Modification 09/06/2013:
# 1. only paired-aligned reads will be reported: with tophat parameter --no-discordant and --no-mixed
# 2. delete the sam files: .sam.all, .sam.minus and .sam.plus 

## Modification:
# 1. the extension of the fastq files will not be fixed in order to run other groups' files

## Modification: v2
# 1. using gencode-v16
# 2. adding spike-ins as reference genome
# 3. merge accepted_hits.bam and unmapped.bam
# 4. novel-juncs only

# usage : submitJob perl /home/CAM/xwei/tgchome/script/encode_pipeline/encode_tophat_pipeline_singleEnd_unstranded.pl $infile_folder $outfile_folder $file_1 $novel_juncs_detection $gender

## --- Proton data ---
# examle: submitJob16 perl /home/CAM/xwei/tgchome/script/encode_pipeline/encode_tophat_pipeline_singleEnd_unstranded.pl /archive/TGCore_User_Data/Proton/evaluation/fastq /archive/TGCore_User_Data/Proton/evaluation/tophat/AKap1_LV08_3 AKap1_LV08_3.IonXpressRNA_003.fastq.gz
##       498328 ("NGS_app")
# examle: submitJob16 perl /home/CAM/xwei/tgchome/script/encode_pipeline/encode_tophat_pipeline_singleEnd_unstranded.pl /archive/TGCore_User_Data/Proton/evaluation/fastq /archive/TGCore_User_Data/Proton/evaluation/tophat/AKap1_LV08_4 AKap1_LV08_4.IonXpressRNA_004.fastq.gz
##       498329 ("NGS_app")






use strict;
use FileHandle;

my $IN_DIR = $ARGV[0]; # the input file dir
my $FILEDIR = $ARGV[1]; # the output dir, e.g. /home/TGC/bgraveley/120511_HWI-EAS299_00035_FC/Project_BGraveley_YeastA_1_0512/Sample_yeast_A3/
my $FASTQ1 = $ARGV[2]; #e.g. 328_AGTTCC_L004_R1_001.fastq.gz  or ENCFF001REI.txt.gz
#my $FASTQ2 = $ARGV[3]; #e.g. 328_AGTTCC_L004_R2_001.fastq.gz  or ENCFF001REJ.txt.gz
my $novel_juncs_detection = $ARGV[4]; #e.g. no-novel-juncs
my $gender = $ARGV[5];  ## e.g. male

my %transcriptome_index_type = ("female" => "/archive/TGCore_User_Data/ENCODE/ENCODE_ANNOTATIONS/gencode.v19.female.spikeins", 
                                "male"   => "/archive/TGCore_User_Data/ENCODE/ENCODE_ANNOTATIONS/gencode.v19.male.spikeins"); 
my %bowtie_index_type = ("female" => "/archive/TGCapps_data/bowtie2_indexes/female.spikeins.hg19 ", 
                         "male"   => "/archive/TGCapps_data/bowtie2_indexes/male.spikeins.hg19 "); 
#my $library_type = "fr-firststrand";
my $library_type = "fr-unstranded"; ## for unstranded
my $chromosizes_file = "/archive/TGCore_User_Data/ENCODE/ENCODE_ANNOTATIONS/chromosizes_spikeIns.txt"; ## add spike ins for the chromosize file

my $fastq1 = $IN_DIR."/".$FASTQ1;
my $fastq2 = $IN_DIR."/".$FASTQ2;

my $output = "";

print "\nstart------\n";
if(! (-d $FILEDIR)){`mkdir $FILEDIR`;}
## 1. Align reads using tophat
my $bowtie2_index = $bowtie_index_type{"female"}; # female index for cell lines K562, GM12878 and Mcf78. male index for cell line HepG2
my $transcriptome_index = $transcriptome_index_type{"female"}; # female index for cell lines K562, GM12878 and Mcf78. male index for cell line HepG2
#if($FILEDIR =~ m/HepG2/ || $FASTQ1 =~ m/HepG2/){
if($gender eq "male"){
    print "using \'male\' reference files\n";
    $bowtie2_index = $bowtie_index_type{"male"};
    $transcriptome_index = $transcriptome_index_type{"male"};
}
###------- special for set4 and set5: use the male genome since the tissues are a mixture of male and females. 
#    $bowtie2_index = $bowtie_index_type{"male"};
#    $transcriptome_index = $transcriptome_index_type{"male"};
###------- remove it when running on cell lines
    
my $no_novel_juncs = "";
if($novel_juncs_detection =~ m /no\-novel\-juncs/){
    $no_novel_juncs = "--no-novel-juncs";
}
print "tophat -p 16 -z0 -o $FILEDIR -a 8 -m 0 --min-intron-length 20 --max-intron-length 1000000 --read-edit-dist 4 --read-mismatches 4 -g 20 $no_novel_juncs --library-type $library_type --transcriptome-index $transcriptome_index $bowtie2_index $fastq1\n";
$output = `tophat -p 16 -z0 -o $FILEDIR -a 8 -m 0 --min-intron-length 20 --max-intron-length 1000000 --read-edit-dist 4 --read-mismatches 4 -g 20 $no_novel_juncs --library-type $library_type --transcriptome-index $transcriptome_index $bowtie2_index $fastq1`;    
print "$output\n";
print "\n";


## 2. Count the tophat aligned reads on accepted_hit.bam file 
my $accepted_file = $FILEDIR . "/accepted_hits.bam";
my $count_file = $accepted_file . ".count";
print "samtools view $accepted_file | wc -l > $count_file\n";
`samtools view $accepted_file | wc -l > $count_file`;
print "\n";


## add more header lines to .bam files
my $unmapped_file = $FILEDIR."/unmapped.bam";
`perl /home/CAM/xwei/tgchome/script/encode_pipeline/tophat_bam_header_revise.pl $accepted_file`;
`perl /home/CAM/xwei/tgchome/script/encode_pipeline/tophat_bam_header_revise.pl $unmapped_file`;

#### index the bam file
my $bam_idx = $accepted_file . ".bai";
print "samtools index $accepted_file $bam_idx\n";
`samtools index $accepted_file $bam_idx`;


## No Need for unstranded samples :: 3. fix the XS:A:+ or XS:A:- tags for tophat version 2.0.8 
my $bam_file = $FILEDIR."/accepted_hits.bam";
my $bam_all = $bam_file;
#my $bam_all = $FILEDIR . "/accepted_hits.all.bam";
#print "perl /home/CAM/xwei/tgchome/script/encode_pipeline/tophat_bam_xsA_tag_fix.pl $bam_file $bam_all\n";
#$output = `perl /home/CAM/xwei/tgchome/script/encode_pipeline/tophat_bam_xsA_tag_fix.pl $bam_file $bam_all`;
#print "$output\n";
###### count the unique mapping reads using tag NH:i:1
my $count_unique_file = $bam_all . ".unique.count";
print "samtools view $bam_all | grep \"NH:i:1\" | wc -l > $count_unique_file\n";
`samtools view $bam_all | grep "NH:i:1" | wc -l > $count_unique_file`;
print "\n";


#### merge the accepted_hits.all.bam and unmapped.bam with the correct +/- tags
my $merged_bam = $FILEDIR."/merged.bam";
print "samtools merge $merged_bam $bam_all $unmapped_file\n";
$output = `samtools merge $merged_bam $bam_all $unmapped_file`;
print "$output\n";


### 4. using Georgi's script (makewigglefromBAM-NH.py) to generate unstranded unique or multi mapped file
##      then convert to bigWig
print "perl /home/CAM/xwei/tgchome/script/encode_pipeline/encode_bam_to_bigwig.pl $bam_file $FILEDIR unstranded unique\n";
$output = `perl /home/CAM/xwei/tgchome/script/encode_pipeline/encode_bam_to_bigwig.pl $bam_file $FILEDIR unstranded unique`;
print "$output\n";
#print "perl /home/CAM/xwei/tgchome/script/encode_pipeline/encode_bam_to_bigwig.pl $bam_file $FILEDIR minus unique\n";
#$output = `perl /home/CAM/xwei/tgchome/script/encode_pipeline/encode_bam_to_bigwig.pl $bam_file $FILEDIR minus unique`;
#print "$output\n";
print "perl /home/CAM/xwei/tgchome/script/encode_pipeline/encode_bam_to_bigwig.pl $bam_file $FILEDIR unstranded multi\n";
$output = `perl /home/CAM/xwei/tgchome/script/encode_pipeline/encode_bam_to_bigwig.pl $bam_file $FILEDIR unstranded multi`;
print "$output\n";
#print "perl /home/CAM/xwei/tgchome/script/encode_pipeline/encode_bam_to_bigwig.pl $bam_file $FILEDIR minus multi\n";
#$output = `perl /home/CAM/xwei/tgchome/script/encode_pipeline/encode_bam_to_bigwig.pl $bam_file $FILEDIR minus multi`;
#print "$output\n";


print "------done\n\n";

####################################
sub get_count{
    my ($count_file) = @_;
    open(COU, $count_file) or die("can't open $count_file"); #INPUT FILE
    my @list = <COU>;
    my $count = pop @list;
    $count =~ s/\n//g;
    my @t = split / /, $count;
    $count = $t[0];
    return $count;
}



