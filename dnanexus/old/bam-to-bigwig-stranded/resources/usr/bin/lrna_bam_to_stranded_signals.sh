#!/bin/bash -e

if [ $# -ne 2 ]; then
    echo "usage v1: lrna_bam_to_stranded_signals.sh <bam_file> <chrom_sizes>"
    echo "Converts BAMs from alignments from stranded libraries to bigwig format. Is independent of DX and encodeD."
    exit -1; 
fi
bam_file=$1      # Bam file.
chrom_sizes=$2   # chrom_sizes file that matches the genome used to create bam_root.

bam_root=${bam_file%.bam}
echo "-- Results will be: '${bam_root}_minusAll.bw', '${bam_root}_minusUniq.bw', '${bam_root}_plusAll.bw', and '${bam_root}_plusUniq.bw'"

echo "-- Make signals..."
set -x
mkdir -p Signal
STAR --runMode inputAlignmentsFromBAM --inputBAMfile $bam_file --outWigType bedGraph \
     --outWigStrand Stranded --outFileNamePrefix ./Signal/ --outWigReferencesPrefix chr
mv Signal/Signal*bg .
set +x

echo "-- Convert bedGraph to bigWigs..."
set -x
bedGraphToBigWig Signal.UniqueMultiple.str1.out.bg $chrom_sizes ${bam_root}_minusAll.bw
bedGraphToBigWig Signal.Unique.str1.out.bg         $chrom_sizes ${bam_root}_minusUniq.bw
bedGraphToBigWig Signal.UniqueMultiple.str2.out.bg $chrom_sizes ${bam_root}_plusAll.bw
bedGraphToBigWig Signal.Unique.str2.out.bg         $chrom_sizes ${bam_root}_plusUniq.bw
set +x

echo "-- The results..."
ls -l ${bam_root}*.bw

