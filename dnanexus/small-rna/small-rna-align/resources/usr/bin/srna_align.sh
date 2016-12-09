#!/bin/bash -e

if [ $# -lt 5 ] ||  [ $# -gt 6 ]; then
    echo "usage v1: srna_align.sh <star_index.tgz> <reads.fq.gz> <library_id> <ncpus> <bam_root> [<clipping_model>]"
    echo "Align single-end short-RNA-seq reads with STAR.  Is independent of DX and encodeD."
    exit -1; 
fi
star_index_tgz=$1  # STAR Index archive.
reads_fq_gz=$2     # gzipped fastq of of single-end reads.
library_id=$3      # Library identifier which will be added to bam header.
ncpus=$4           # Number of cpus available.
bam_root="$5_srna_star" # root name for output bam (e.g. "out_bam" will create "out_bam_srna_star.bam")
clipping_model="ENCODE3"
if [ $# -eq 6 ]; then
    clipping_model=$6  # "A_Tailing_No_Barcode", "A_Tailing_N3", "A_Tailing_N4"
fi

echo "-- Alignments file will be: '${bam_root}.bam'"

echo "-- Extracting star index archive..."
tar zxvf $star_index_tgz
# unzips into "out/"

echo "-- Set up headers..."
set -x
libraryComment="@CO\tLIBID:${library_id}"
echo -e ${libraryComment} > COfile.txt
cat out/*_bamCommentLines.txt >> COfile.txt
echo `cat COfile.txt`
set +x

clip_params="--clip3pAdapterSeq TGGAATTCTC --clip3pAdapterMMp 0.1"
if [ $clipping_model == "A_Tailing_No_Barcode" ]; then
    clip_params="--clip3pAdapterSeq AAAAA --clip3pAdapterMMp 0.0 --clip5pNbases 0"
elif [ $clipping_model == "A_Tailing_N3" ]; then
    clip_params="--clip3pAdapterSeq AAAAA --clip3pAdapterMMp 0.0 --clip5pNbases 5"
elif [ $clipping_model == "A_Tailing_N4" ]; then
    clip_params="--clip3pAdapterSeq AAAAA --clip3pAdapterMMp 0.0 --clip5pNbases 6"
elif [ $clipping_model != "ENCODE3" ]; then
    echo "-- WARNING: Unknown clipping model '$clipping_model'"
    clipping_model="ENCODE3"
fi
echo "-- Using clipping model '$clipping_model'"


echo "-- Map reads..."
set -x
STAR --genomeDir out --readFilesIn $reads_fq_gz --readFilesCommand zcat                     \
    --runThreadN $ncpus --outFilterMultimapNmax 20 --alignIntronMax 1                       \
    $clip_params --outFilterMismatchNoverLmax 0.03                                          \
    --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 16  \
    --outSAMheaderCommentFile COfile.txt --outSAMheaderHD @HD VN:1.4 SO:coordinate          \
    --genomeLoad NoSharedMemory --outSAMunmapped Within --outSAMtype BAM SortedByCoordinate \
    --quantMode GeneCounts --alignSJDBoverhangMin 1000 --limitBAMsortRAM 60000000000
        
mv Aligned.sortedByCoord.out.bam ${bam_root}.bam
mv ReadsPerGene.out.tab ${bam_root}_quant.tsv
mv Log.final.out ${bam_root}_Log.final.out
set +x

echo "-- Collect bam flagstats..."
set -x
samtools flagstat ${bam_root}.bam > ${bam_root}_flagstat.txt
set +x

echo "-- The results..."
ls -l ${bam_root}*

