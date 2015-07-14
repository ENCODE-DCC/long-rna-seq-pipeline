#!/usr/bin/perl
#
# tophat_bam_header_revise.pl
# 
# xwei  7/6/2013

## add more information to the header of .bam files

## example:   perl /home/CAM/xwei/tgchome/script/encode_pipeline/tophat_bam_header_revise.pl /home/CAM/xwei/tgchome/tophatTest_15/YEO_8-16/unmapped.bam



use strict;
use FileHandle;

my $BAM= $ARGV[0];

my $add_header_file = "/home/ubuntu/data/long-rna-seq-pipeline/DAC/additional_bam_header.txt";

my $header_file = $BAM . ".header";
my $new_bam = $BAM;
$new_bam =~ s/bam/newheader\.bam/;

`samtools view -H $BAM > $header_file`;
`cat $add_header_file >> $header_file`;
`samtools reheader $header_file $BAM > $new_bam`;
`rm $BAM`;
`mv $new_bam $BAM`;
`rm $header_file`;






