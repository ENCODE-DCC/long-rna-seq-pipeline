#!/bin/bash -e

if [ $# -ne 5 ]; then
    echo "usage v1: lrna-align-star-se.sh <star_index.tgz> <reads.fq.gz> <library_id> <ncpus> <bam_root>"
    echo "Align single-end reads with STAR.  Is independent of DX and encodeD."
    exit -1; 
fi
star_index_tgz=$1  # STAR Index archive.
reads_fq_gz=$2     # gzipped fastq of of single-end reads.
library_id=$3      # Library identifier which will be added to bam header.
ncpus=$4            # Number of cpus available.
bam_root="$5_star" # root name for output bam (e.g. "out_bam" will create "out_bam_star_genome.bam" and "out_bam_star_anno.bam") 

echo "-- Alignments file will be: '${bam_root}_genome.bam' and '${bam_root}_anno.bam'"

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
STAR --genomeDir out --readFilesIn $reads_fq_gz                                 \
    --readFilesCommand zcat --runThreadN $ncpus --genomeLoad NoSharedMemory      \
    --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1    \
    --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04              \
    --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000         \
    --outSAMheaderCommentFile COfile.txt --outSAMheaderHD @HD VN:1.4 SO:coordinate   \
    --outSAMunmapped Within --outFilterType BySJout --outSAMattributes NH HI AS NM MD \
    --outSAMstrandField intronMotif --outSAMtype BAM SortedByCoordinate                \
    --quantMode TranscriptomeSAM --sjdbScore 1 --limitBAMsortRAM 60000000000

mv Aligned.sortedByCoord.out.bam ${bam_root}_genome.bam
mv Log.final.out ${bam_root}_Log.final.out
set +x
ls -l ${bam_root}_genome.bam

echo "-- Sorting annotation bam..."
set -x
cat <( samtools view -H Aligned.toTranscriptome.out.bam ) \
    <( samtools view -@ $ncpus Aligned.toTranscriptome.out.bam | sort -S 60G -T ./ ) | \
    samtools view -@ $ncpus -bS - > ${bam_root}_anno.bam
set +x
ls -l ${bam_root}_anno.bam

echo "-- Collect bam flagstats..."
set -x
samtools flagstat ${bam_root}_genome.bam > ${bam_root}_genome_flagstat.txt
samtools flagstat ${bam_root}_anno.bam > ${bam_root}_anno_flagstat.txt
set +x

echo "-- The results..."
ls -l ${bam_root}*

