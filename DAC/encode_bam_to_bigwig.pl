#!/usr/bin/perl
#
# encode_bam_to_bigwig.pl
# 
# xwei  09/24/2013

### modification 9/24/2013
# 1. can generate unstranded wig files

### using Georgi's script (makewigglefromBAM-NH.py) to split accepted_hits.bam to plus or minus strand. and generate unique or multi mapped files

### usage: /home/CAM/xwei/tgchome/script/encode_pipeline/encode_bam_to_bigwig.pl <accepted_hits.bam> <output folder> <plus | minus> <unique | multi>

### example : /home/CAM/xwei/tgchome/script/encode_pipeline/encode_bam_to_bigwig.pl /archive/TGCore_User_Data/ENCODE/ENCODE_DATA/set1/bigwig/test_make3/accepted_hits.bam /archive/TGCore_User_Data/ENCODE/ENCODE_DATA/set1/bigwig/test_make3 plus unique
### example : /home/CAM/xwei/tgchome/script/encode_pipeline/encode_bam_to_bigwig.pl /archive/TGCore_User_Data/ENCODE/ENCODE_DATA/set1/bigwig/test_make3/accepted_hits.bam /archive/TGCore_User_Data/ENCODE/ENCODE_DATA/set1/bigwig/test_make3 minus multi

### submitJob1 /home/CAM/xwei/tgchome/script/encode_pipeline/encode_bam_to_bigwig.pl /archive/TGCore_User_Data/ENCODE/ENCODE_DATA/set1/tophat/novel_juncs/GM12878-Br-1/accepted_hits.bam /archive/TGCore_User_Data/ENCODE/ENCODE_DATA/set1/tophat/novel_juncs/GM12878-Br-1 minus multi     
### submitJob1 /home/CAM/xwei/tgchome/script/encode_pipeline/encode_bam_to_bigwig.pl /archive/TGCore_User_Data/ENCODE/ENCODE_DATA/set1/tophat/novel_juncs/GM12878-Br-1/accepted_hits.bam /archive/TGCore_User_Data/ENCODE/ENCODE_DATA/set1/tophat/novel_juncs/GM12878-Br-1 plus multi     
### submitJob1 /home/CAM/xwei/tgchome/script/encode_pipeline/encode_bam_to_bigwig.pl /archive/TGCore_User_Data/ENCODE/ENCODE_DATA/set1/tophat/novel_juncs/GM12878-Br-1/accepted_hits.bam /archive/TGCore_User_Data/ENCODE/ENCODE_DATA/set1/tophat/novel_juncs/GM12878-Br-1 plus unique     
### submitJob1 /home/CAM/xwei/tgchome/script/encode_pipeline/encode_bam_to_bigwig.pl /archive/TGCore_User_Data/ENCODE/ENCODE_DATA/set1/tophat/novel_juncs/GM12878-Br-1/accepted_hits.bam /archive/TGCore_User_Data/ENCODE/ENCODE_DATA/set1/tophat/novel_juncs/GM12878-Br-1 minus unique     

### submitJob1 /home/CAM/xwei/tgchome/script/encode_pipeline/encode_bam_to_bigwig.pl /archive/TGCore_User_Data/ENCODE/ENCODE_DATA/set1/tophat/novel_juncs/GM12878-Br-2/accepted_hits.bam /archive/TGCore_User_Data/ENCODE/ENCODE_DATA/set1/tophat/novel_juncs/GM12878-Br-2 minus multi     
### submitJob1 /home/CAM/xwei/tgchome/script/encode_pipeline/encode_bam_to_bigwig.pl /archive/TGCore_User_Data/ENCODE/ENCODE_DATA/set1/tophat/novel_juncs/GM12878-Br-2/accepted_hits.bam /archive/TGCore_User_Data/ENCODE/ENCODE_DATA/set1/tophat/novel_juncs/GM12878-Br-2 plus multi     
### submitJob1 /home/CAM/xwei/tgchome/script/encode_pipeline/encode_bam_to_bigwig.pl /archive/TGCore_User_Data/ENCODE/ENCODE_DATA/set1/tophat/novel_juncs/GM12878-Br-2/accepted_hits.bam /archive/TGCore_User_Data/ENCODE/ENCODE_DATA/set1/tophat/novel_juncs/GM12878-Br-2 plus unique     
### submitJob1 /home/CAM/xwei/tgchome/script/encode_pipeline/encode_bam_to_bigwig.pl /archive/TGCore_User_Data/ENCODE/ENCODE_DATA/set1/tophat/novel_juncs/GM12878-Br-2/accepted_hits.bam /archive/TGCore_User_Data/ENCODE/ENCODE_DATA/set1/tophat/novel_juncs/GM12878-Br-2 minus unique     

### submitJob1 /home/CAM/xwei/tgchome/script/encode_pipeline/encode_bam_to_bigwig.pl /archive/TGCore_User_Data/ENCODE/ENCODE_DATA_OTHERS/gi_set1/tophat/novel_juncs/GM12878_Long_PolyA_dUTP_2x100_ENCLB038ZZZ_Ging/accepted_hits.bam /archive/TGCore_User_Data/ENCODE/ENCODE_DATA_OTHERS/gi_set1/tophat/novel_juncs/GM12878_Long_PolyA_dUTP_2x100_ENCLB038ZZZ_Ging minus multi     
#32959.tgccluster.tgc.cluster
### submitJob1 /home/CAM/xwei/tgchome/script/encode_pipeline/encode_bam_to_bigwig.pl /archive/TGCore_User_Data/ENCODE/ENCODE_DATA_OTHERS/gi_set1/tophat/novel_juncs/GM12878_Long_PolyA_dUTP_2x100_ENCLB038ZZZ_Ging/accepted_hits.bam /archive/TGCore_User_Data/ENCODE/ENCODE_DATA_OTHERS/gi_set1/tophat/novel_juncs/GM12878_Long_PolyA_dUTP_2x100_ENCLB038ZZZ_Ging plus multi     
#32958.tgccluster.tgc.cluster
### submitJob1 /home/CAM/xwei/tgchome/script/encode_pipeline/encode_bam_to_bigwig.pl /archive/TGCore_User_Data/ENCODE/ENCODE_DATA_OTHERS/gi_set1/tophat/novel_juncs/GM12878_Long_PolyA_dUTP_2x100_ENCLB038ZZZ_Ging/accepted_hits.bam /archive/TGCore_User_Data/ENCODE/ENCODE_DATA_OTHERS/gi_set1/tophat/novel_juncs/GM12878_Long_PolyA_dUTP_2x100_ENCLB038ZZZ_Ging plus unique     
#32956.tgccluster.tgc.cluster
### submitJob1 /home/CAM/xwei/tgchome/script/encode_pipeline/encode_bam_to_bigwig.pl /archive/TGCore_User_Data/ENCODE/ENCODE_DATA_OTHERS/gi_set1/tophat/novel_juncs/GM12878_Long_PolyA_dUTP_2x100_ENCLB038ZZZ_Ging/accepted_hits.bam /archive/TGCore_User_Data/ENCODE/ENCODE_DATA_OTHERS/gi_set1/tophat/novel_juncs/GM12878_Long_PolyA_dUTP_2x100_ENCLB038ZZZ_Ging minus unique     
#32957.tgccluster.tgc.cluster


use strict;
use FileHandle;

my $BAM = $ARGV[0];
my $OUTDIR = $ARGV[1];
my $strand_type = $ARGV[2]; ## puls or minus
my $map_type = $ARGV[3]; ## unique or multi

my $makewigglefromBAM = "/home/CAM/xwei/tgchome/script/Marinov_Georgi/V5_09192013/makewigglefromBAM-NH.py"; #   V4_09072013  V5_09192013 -- plus tracks are the same as minus tracks
my $chromosizes = "/archive/TGCore_User_Data/ENCODE/ENCODE_ANNOTATIONS/chromosizes.txt"; ## no spikeIns


print "\n-- start --\n";
#### 1. index the bam fle first: samtools index test_sorted.bam test_sorted.bai
##example: submitJob "samtools index /archive/TGCore_User_Data/ENCODE/ENCODE_DATA/set1/bigwig/test_make2/accepted_hits.bam /archive/TGCore_User_Data/ENCODE/ENCODE_DATA/set1/bigwig/test_make2/accepted_hits.bam.bai"   
#my $bam_idx = $BAM . ".bai";
#print "samtools index $BAM $bam_idx\n";
##`samtools index $BAM $bam_idx`;   

### 2. convert files
my $common_option = "-RPM -notitle -fragments second-read-strand";
my $strand = "-stranded";
my $out_file = $OUTDIR;
if($strand_type eq "plus"){
    $strand .= " +";
    $out_file .= "/plus.";
}elsif($strand_type eq "minus"){
    $strand .= " -";
    $out_file .= "/minus.";
}elsif($strand_type eq "unstranded"){
    $strand = "";
    $out_file .= "/";
}
my $mapping = "";
if($map_type eq "unique"){
    $mapping = "-nomulti";
    $out_file .= "unique";
}else{
    $out_file .= "multi";
}
$out_file .= ".scaled.2";

my $wig = $out_file . ".wig";
my $rm_wig = $out_file . ".rm.wig";
my $bw = $out_file . ".bw";

# (1) bam to wig
print "python $makewigglefromBAM --- $BAM $chromosizes $wig $strand $mapping $common_option \n";
`python $makewigglefromBAM --- $BAM $chromosizes $wig $strand $mapping $common_option`;
# (2) remove "-"sign from minus files
if($strand_type eq "plus"){
    $rm_wig = $wig;
}elsif($strand_type eq "minus"){
    print "perl -pe 's/\t\-/\t/g' < $wig > $rm_wig\n";
    `perl -pe 's/\t\-/\t/g' < $wig > $rm_wig`;
}elsif($strand_type eq "unstranded"){
    $rm_wig = $wig;
}
# (3) wigToBigWig
print "wigToBigWig $rm_wig $chromosizes $bw\n";
`wigToBigWig $rm_wig $chromosizes $bw`;

print "\n\n";






########################################################
#### 2. run the script
## example: submitJob "python /home/CAM/xwei/tgchome/script/Marinov_Georgi/V4_09072013/makewigglefromBAM-NH.py --- /archive/TGCore_User_Data/ENCODE/ENCODE_DATA/set1/bigwig/test_make2/accepted_hits.bam /archive/TGCore_User_Data/ENCODE/ENCODE_DATA/set1/bigwig/test_make2/chromosizes_spikeIns.txt /archive/TGCore_User_Data/ENCODE/ENCODE_DATA/set1/bigwig/test_make2/v4.plus.unique.scaled.2.wig -stranded + -nomulti -RPM -notitle -fragments second-read-strand "   
#my $plus_uni_wig = $OUTDIR . "/plus.unique.scaled.2.wig";
#print "python $makewigglefromBAM --- $BAM $chromosizes $plus_uni_wig -stranded + -nomulti -RPM -notitle -fragments second-read-strand \n";   
#my $minus_uni_wig = $OUTDIR . "/minus.unique.scaled.2.wig";
#print "python $makewigglefromBAM --- $BAM $chromosizes $minus_uni_wig -stranded - -nomulti -RPM -notitle -fragments second-read-strand \n";   
#my $plus_mul_wig = $OUTDIR . "/plus.multi.scaled.2.wig";
#print "python $makewigglefromBAM --- $BAM $chromosizes $plus_mul_wig -stranded + -RPM -notitle -fragments second-read-strand \n";   
#my $minus_mul_wig = $OUTDIR . "/minus.multi.scaled.2.wig";
#print "python $makewigglefromBAM --- $BAM $chromosizes $minus_mul_wig -stranded - -RPM -notitle -fragments second-read-strand \n"; 
#
#### 3. remove "-"sign from minus files
## example: submitJob1 "perl -pe 's/\t\-/\t/g' < /archive/TGCore_User_Data/ENCODE/ENCODE_DATA/set1/bigwig/test_make2/v4.minus.unique.scaled.2.wig > /archive/TGCore_User_Data/ENCODE/ENCODE_DATA/set1/bigwig/test_make2/v4.minus.unique.scaled.2.rm.wig"
#my $minus_uni_rm_wig = $minus_uni_wig;
#$minus_uni_rm_wig =~ s /\.wig/\.rm\.wig/;
#print "perl -pe 's/\t\-/\t/g' < $minus_uni_wig > $minus_uni_rm_wig\n";
#my $minus_mul_rm_wig = $minus_mul_wig;
#$minus_mul_rm_wig =~ s /\.wig/\.rm\.wig/;
#print "perl -pe 's/\t\-/\t/g' < $minus_mul_wig > $minus_mul_rm_wig\n";
#
#### 4. wigToBigWig
## example: submitJob1 "wigToBigWig /archive/TGCore_User_Data/ENCODE/ENCODE_DATA/set1/bigwig/test_make2/v4.minus.unique.scaled.2.rm.wig /archive/TGCore_User_Data/ENCODE/ENCODE_ANNOTATIONS/chromosizes_spikeIns.txt  /archive/TGCore_User_Data/ENCODE/ENCODE_DATA/set1/bigwig/test_make2/v4.minus.unique.scaled.2.bw"
#my $plus_uni_bw = $plus_uni_wig;
#$plus_uni_bw =~ s /\.wig/\.bw/;
#print "wigToBigWig $plus_uni_wig $chromosizes $plus_uni_bw\n";
#my $minus_uni_bw = $minus_uni_wig;
#$minus_uni_bw =~ s /\.wig/\.bw/;
#print "wigToBigWig $minus_uni_rm_wig $chromosizes $minus_uni_bw\n";
#my $plus_mul_bw = $plus_mul_wig;
#$plus_mul_bw =~ s /\.wig/\.bw/;
#print "wigToBigWig $plus_mul_wig $chromosizes $plus_mul_bw\n";
#my $minus_mul_bw = $minus_mul_wig;
#$minus_mul_bw =~ s /\.wig/\.bw/;
#print "wigToBigWig $minus_mul_rm_wig $chromosizes $minus_mul_bw\n";






