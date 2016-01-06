#!/bin/bash -e

if [ $# -ne 2 ]; then
    echo "usage v1: lrna_bam_to_unstranded_signals.sh <bam_file> <chrom_sizes>"
    echo "Converts BAMs from alignments from unstranded libraries to bigwig format. Is independent of DX and encodeD."
    exit -1; 
fi
bam_file=$1      # Bam file.
chrom_sizes=$2   # chrom_sizes file that matches the genome used to create bam_root.

bam_root=${bam_file%.bam}
echo "-- Results will be: '${bam_root}_all.bw' and '${bam_root}_uniq.bw'"

echo "-- Make signals..."
set -x
mkdir -p Signal
STAR --runMode inputAlignmentsFromBAM --inputBAMfile $bam_file --outWigType bedGraph \
     --outWigStrand Unstranded --outFileNamePrefix ./Signal/ --outWigReferencesPrefix chr
mv Signal/Signal*bg .
set +x
echo `ls -l`

echo "-- Convert bedGraph to bigWigs..."
set -x
bedGraphToBigWig Signal.UniqueMultiple.str1.out.bg chromSizes.txt ${bam_root}_all.bw
bedGraphToBigWig Signal.Unique.str1.out.bg         chromSizes.txt ${bam_root}_uniq.bw
set +x

echo "-- The results..."
ls -l ${bam_root}*.bw

