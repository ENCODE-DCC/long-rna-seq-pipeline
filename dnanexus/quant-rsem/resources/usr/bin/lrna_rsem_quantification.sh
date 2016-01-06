#!/bin/bash -e

if [ $# -ne 5 ]; then
    echo "usage v1: lrna_rsem_quantification.sh <rsem_index.tgz> <anno_bam> <paired_end> <rnd_seed> <ncpus>"
    echo "Align single-end reads with STAR.  Is independent of DX and encodeD."
    exit -1; 
fi
rsem_index_tgz=$1  # RSEM Index archive.
anno_bam=$2        # STAR alignment to annotation
paired_end=$3      # "true" if alignment was on paired-end data.
rnd_seed=$4        # Random seed.  ENCODE has been using 12345
ncpus=$5            # Number of cpus available.

bam_root=${anno_bam%.bam}
echo "-- Qunatification results will be: '${bam_root}_rsem.genes.results' and '${bam_root}_rsem.isoforms.results'"

echo "-- Extracting star index archive..."
tar zxvf $rsem_index_tgz
# should be 'out/rsem'

grp=`ls out/*.grp`
index_prefix=${grp%.grp}
echo "-- Found index_prefix: '$index_prefix'"

extraFlags=""
if [ "$paired_end" == "true" ]; then
    echo '-- Running for paired-end, stranded'
    extraFlags="--paired-end --forward-prob 0"
else
    echo '-- Running for unpaired, unstranded'
fi

#if [ "$stranded" == "true" ]
#then
#    echo '* Using stranded flag'
#    extraFlags=${extraFlags}"--forward-prob 0"
#fi

echo "-- Quantify with extra flags: [${extraFlags}]..."
set -x
rsem-calculate-expression --bam --estimate-rspd --calc-ci --seed ${rnd_seed} -p $ncpus \
    --no-bam-output --ci-memory 30000 ${extraFlags} $anno_bam ${index_prefix} ${bam_root}_rsem
set +x

echo "-- The results..."
ls -l ${bam_root}*.results

