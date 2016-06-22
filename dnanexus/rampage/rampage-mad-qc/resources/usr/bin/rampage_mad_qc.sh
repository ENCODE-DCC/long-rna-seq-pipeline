#!/bin/bash -e

if [ $# -lt 3 ] || [ $# -gt 4 ]; then
    echo "usage v1: rampage_mad_qc.sh <quants_file_a> <quants_file_b> <out_root_name> [MAD_R_script]"
    echo "Generages MAD QC scoring of 2 replicate quantification files and produces: *_plot.png/.json. Is independent of DX and encodeD."
    exit -1; 
fi
quants_file_a=$1     # One replicate quantification file 
quants_file_b=$2     # Another replicate quantification file
mad_root_name=$3     # Root name for output ((e.g. 'mad' will result in mad_plot.png, and mad.json) 
MAD_R_script='/usr/bin/MAD.R'
if [ $# -eq 4 ]; then
    MAD_R_script=$4          # Path to R_script (e.g. '/usr/bin/MAD.R')
fi 


echo "-- Running '$MAD_R_script'..."
Rscript $MAD_R_script $quants_file_a $quants_file_b > ${mad_root_name}.json
mv MAplot.png ${mad_root_name}_plot.png

echo "-- The results..."
ls -l ${mad_root_name}*

