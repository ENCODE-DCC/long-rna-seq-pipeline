#!/bin/bash -e

if [ $# -ne 7 ]; then
    echo "usage v1: rampage_peaks.sh <rampage_bam> <control_bam> <annotation_gtf_gz> <chrom_sizes> <'rampage'|'cage'> <ncpus> <peaks_root_name>"
    echo "Generages Rampage Peaks from marked bam, producing files: *.bed.gz/.bb/.gtf.gz/.tsv. Is independent of DX and encodeD."
    exit -1; 
fi
rampage_bam=$1       # Rampage dumplicates marked bam file.
control_bam=$2       # long-RNA-seq control bam file.
annotation_gtf_gz=$3 # Gzipped annotations file in gtf format.
chrom_sizes=$4       # chrom_sizes file that matches the genome used to create bam_root.
assay_type=$5        # Is this cage or rampage data?
ncpus=$6             # Number of cpus available.
peaks_root=$7        # Root name of peaks result file (e.g. 'peaks' will result in peaks.bed.gz, peaks.bb, peaks.gtf.gz and peaks.tsv) 

if [ "$assay_type" == "rampage" ] || [ "$assay_type" == "cage" ] ; then
    echo '-- Running for $assay_type reads'
    reads_flag="--${assay_type}-reads"
else
    echo '-- Unknown assay type: $assay_type' 
    exit -1; 
fi

bam_root=${rampage_bam%.bam}
echo "-- Peaks files will be: '${peaks_root}.*'"

annotation_file=$annotation_gtf_gz
if [[ "$annotation_file" == *.gz ]]; then
    echo "-- Uncompressing ${annotation_file} file..."
    set -x
    gunzip $annotation_file
    set +x
    annotation_file=${annotation_file%.gz}
fi 

echo "-- Indexing bams if necessary..."
if [ ! -e "${rampage_bam}.bai" ]; then
    set -x
    samtools index $rampage_bam
    set +x
fi 
if [ ! -e "${control_bam}.bai" ]; then
    set -x
    samtools index $control_bam 
    set +x
fi 

echo "-- Calling peaks..."
set -x
call_peaks $reads_flag $rampage_bam --rnaseq-reads $control_bam --threads $ncpus \
           --reference $annotation_file --exp-filter-fraction 0.05 --trim-fraction 0.01 \
           --ucsc --outfname ${peaks_root}.gff --outfname-type gff \
           --bed-peaks-ofname ${peaks_root}.bed \
           --annotation-quantifications-ofname ${peaks_root}_quant.tsv
set +x
 
echo "-- Removing 'chrphiX' from bed..."
# Because validateFiles fails on too long of bed lines, we now trim column 11 to first 512 bytes
set -x
grep -v "^track" ${peaks_root}.bed | sort -k1,1 -k2,2n > peaks.bed
cat peaks.bed | awk -F '\t' '{printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%.512s\n",$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11}' \
                > ${peaks_root}.bed
set +x

echo "-- Converting bed to bigBed..."
set -x
grep "^chr" ${peaks_root}.bed | sort -k1,1 -k2,2n > peaks_polished.bed
bedToBigBed peaks_polished.bed -type=bed6+ -as=/usr/bin/tss_peak.as $chrom_sizes ${peaks_root}.bb
set +x
 
echo "-- Compressing bed and gff..."
set -x
pigz ${peaks_root}.bed
pigz ${peaks_root}.gff
set +x

echo "-- The results..."
ls -l ${peaks_root}*

