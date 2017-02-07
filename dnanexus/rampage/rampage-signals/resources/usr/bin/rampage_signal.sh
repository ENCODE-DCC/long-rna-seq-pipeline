#!/bin/bash -e

if [ $# -lt 3 ] ||  [ $# -gt 4 ]; then
    echo "usage v1: rampage_signal.sh <bam_file> <chrom_sizes> <signal_root_name> [<stranded:T/F>]"
    echo "Converts BAMs from rampage alignments to signal files: *._rampage_5p.bw. Is independent of DX and encodeD."
    exit -1; 
fi
bam_file=$1      # Bam file.
chrom_sizes=$2   # chrom_sizes file that matches the genome used to create bam_root.
signal_root=$3   # Root name of signals result files (e.g. 'signal' will result in signal_minusAll.bw, signal_minusUniq.bw, signal_plusAll.bw and signal_plusUniq.bw)
stranded=:"F"
if  [ $# -eq 4 ]; then 
    stranded=${4:0:1}  # Strand specific library.  Only need first letter to be T or Y otherwise default unstranded
fi

bam_root=${bam_file%.bam}

echo "-- Signal files will be: '${signal_root}_*.bw'"

# force uppercase in compare
if [ "${stranded^^}" == "T" ] || [ "${stranded^^}" == "Y" ] || [ "${stranded}" == "1" ]; then

    echo "-- Make stranded signals..."
    set -x
    STAR --runMode inputAlignmentsFromBAM --inputBAMfile $bam_file --outWigType bedGraph read1_5p \
        --outWigStrand Stranded --outFileNamePrefix read1_5p. --outWigReferencesPrefix chr
    set +x

    echo "-- Convert stranded bedGraph to bigWigs..."
    set -x
    bedGraphToBigWig read1_5p.Signal.UniqueMultiple.str2.out.bg $chrom_sizes ${signal_root}_minusAll.bw
    bedGraphToBigWig read1_5p.Signal.Unique.str2.out.bg         $chrom_sizes ${signal_root}_minusUniq.bw
    bedGraphToBigWig read1_5p.Signal.UniqueMultiple.str1.out.bg $chrom_sizes ${signal_root}_plusAll.bw
    bedGraphToBigWig read1_5p.Signal.Unique.str1.out.bg         $chrom_sizes ${signal_root}_plusUniq.bw
    set +x

else

    echo "-- Make unstranded signals..."
    set -x
    STAR --runMode inputAlignmentsFromBAM --inputBAMfile $bam_file --outWigType bedGraph read1_5p \
        --outWigStrand Unstranded --outFileNamePrefix read1_5p. --outWigReferencesPrefix chr
    set +x
    echo `ls -l`

    echo "-- Convert unstranded bedGraph to bigWigs..."
    set -x
    bedGraphToBigWig read1_5p.Signal.UniqueMultiple.str1.out.bg $chrom_sizes ${signal_root}_all.bw
    bedGraphToBigWig read1_5p.Signal.Unique.str1.out.bg         $chrom_sizes ${signal_root}_uniq.bw
    set +x

fi

echo "-- The results..."
ls -l ${signal_root}_*.bw

