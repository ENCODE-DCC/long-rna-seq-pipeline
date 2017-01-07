#!/bin/bash -e

if [ $# -lt 2 ] ||  [ $# -gt 3 ]; then
    echo "usage v1: lrna_bam_to_stranded_signals.sh <bam_file> <chrom_sizes> [<stranded:T/F>]"
    echo "Converts BAMs from alignments from stranded or unstranded libraries to bigwig format. Is independent of DX and encodeD."
    exit -1; 
fi
bam_file=$1         # Bam file.
chrom_sizes=$2      # chrom_sizes file that matches the genome used to create bam_root.
stranded=:"F"
if  [ $# -eq 3 ]; then 
    stranded=${3:0:1}  # Strand specific library.  Only need first letter to be T or Y otherwise default unstranded
fi

bam_root=${bam_file%.bam}
echo "-- Results will be: '${bam_root}_*.bw'"

# force uppercase in compare
if [ "${stranded^^}" == "T" ] || [ "${stranded^^}" == "Y" ] || [ "${stranded}" == "1" ]; then

    echo "-- Make stranded signals..."
    set -x
    mkdir -p Signal
    STAR --runMode inputAlignmentsFromBAM --inputBAMfile $bam_file --outWigType bedGraph \
        --outWigStrand Stranded --outFileNamePrefix ./Signal/ --outWigReferencesPrefix chr
    mv Signal/Signal*bg .
    set +x

    echo "-- Convert stranded bedGraph to bigWigs..."
    set -x
    bedGraphToBigWig Signal.UniqueMultiple.str1.out.bg $chrom_sizes ${bam_root}_minusAll.bw
    bedGraphToBigWig Signal.Unique.str1.out.bg         $chrom_sizes ${bam_root}_minusUniq.bw
    bedGraphToBigWig Signal.UniqueMultiple.str2.out.bg $chrom_sizes ${bam_root}_plusAll.bw
    bedGraphToBigWig Signal.Unique.str2.out.bg         $chrom_sizes ${bam_root}_plusUniq.bw
    set +x

else

    echo "-- Make unstranded signals..."
    set -x
    mkdir -p Signal
    STAR --runMode inputAlignmentsFromBAM --inputBAMfile $bam_file --outWigType bedGraph \
        --outWigStrand Unstranded --outFileNamePrefix ./Signal/ --outWigReferencesPrefix chr
    mv Signal/Signal*bg .
    set +x
    echo `ls -l`

    echo "-- Convert unstranded bedGraph to bigWigs..."
    set -x
    bedGraphToBigWig Signal.UniqueMultiple.str1.out.bg $chrom_sizes ${bam_root}_all.bw
    bedGraphToBigWig Signal.Unique.str1.out.bg         $chrom_sizes ${bam_root}_uniq.bw
    set +x

fi

echo "-- The results..."
ls -l ${bam_root}*.bw

