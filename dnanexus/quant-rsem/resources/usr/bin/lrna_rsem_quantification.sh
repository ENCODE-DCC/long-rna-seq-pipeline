#!/bin/bash -e

if [ $# -ne 6 ]; then
    echo "usage v1: lrna_rsem_quantification.sh <rsem_index.tgz> <anno_bam> <paired_end> <read_strand> <rnd_seed> <ncpus>"
    echo "Align single-end reads with STAR.  Is independent of DX and encodeD."
    exit -1; 
fi
rsem_index_tgz=$1  # RSEM Index archive.
anno_bam=$2        # STAR alignment to annotation
paired_end=$3      # "true" if alignment was on paired-end data.
read_strand=$4     # strandedness of read (forward, reverse, unstranded)
rnd_seed=$5        # Random seed.  ENCODE has been using 12345
ncpus=$6            # Number of cpus available.


bam_root=${anno_bam%.bam}
echo "-- Qunatification results will be: '${bam_root}_rsem.genes.results' and '${bam_root}_rsem.isoforms.results'"

echo "-- Extracting star index archive..."
tar zxvf $rsem_index_tgz
# should be 'out/rsem'

grp=`ls out/*.grp`
index_prefix=${grp%.grp}
echo "-- Found index_prefix: '$index_prefix'"

extra_flags=""
msg=""
if [ "$paired_end" == "true" ]; then
    msg="paired-end"
    extra_flags="--paired-end"
else
    msg="single-end"
fi
if [ "$read_strand" == "forward" ] || [ "$read_strand" == "+" ] || [ "$read_strand" == "ScriptSeq" ]; then
    extra_flags="$extra_flags --forward-prob 1"
    msg="$msg forward strand"
elif [ "$read_strand" == "reverse" ] || [ "$read_strand" == "-" ] || [ "$read_strand" == "TruSeq" ]; then
    extra_flags="$extra_flags --forward-prob 0"
    msg="$msg reverse strand"
else
    extra_flags="$extra_flags --forward-prob 0.5"
    msg="$msg unstranded"
fi
echo "-- Running as $msg"

echo "-- Quantify with extra flags: [${extra_flags}]..."
set -x
rsem-calculate-expression --bam --estimate-rspd --calc-ci --seed ${rnd_seed} -p $ncpus \
    --no-bam-output --ci-memory 30000 ${extra_flags} $anno_bam ${index_prefix} ${bam_root}_rsem
set +x

echo "-- The results..."
ls -l ${bam_root}*.results

