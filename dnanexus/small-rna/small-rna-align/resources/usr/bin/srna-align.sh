#!/bin/bash -e

if [ $# -ne 5 ]; then
    echo "usage v1: srna-align.sh <star_index.tgz> <reads.fq.gz> <library_id> <ncpus> <bam_root>"
    echo "Align single-end short-RNA-seq reads with STAR.  Is independent of DX and encodeD."
    exit -1; 
fi
star_index_tgz=$1  # STAR Index archive.
reads_fq_gz=$2     # gzipped fastq of of single-end reads.
library_id=$3      # Library identifier which will be added to bam header.
ncpus=$4           # Number of cpus available.
bam_root="$5_srna_star" # root name for output bam (e.g. "out_bam" will create "out_bam_srna_star.bam") 

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

echo "-- Map reads..."
set -x
STAR --genomeDir out --readFilesIn $reads_fq_gz --readFilesCommand zcat                     \
    --runThreadN $ncpus --outFilterMultimapNmax 20 --alignIntronMax 1                       \
    --clip3pAdapterSeq TGGAATTCTC --clip3pAdapterMMp 0.1 --outFilterMismatchNoverLmax 0.03  \
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

