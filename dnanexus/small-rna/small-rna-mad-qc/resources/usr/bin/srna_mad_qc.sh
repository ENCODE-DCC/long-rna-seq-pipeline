#!/bin/bash -e

if [ $# -ne 4 ]; then
    echo "usage v1: srna-mad-qc.sh <annotation_gtf_gz> <quants_a_tsv> <quants_a_tsv> <out_root>"
    echo "Calculates Mean Absolute Deviation and other stats on a pair of quantifications. Is independent of DX and encodeD."
    echo "Expects extract_gene_ids.awk, sum_srna_expression.awk and MAD.R in current directory."
    exit -1; 
fi
annotation_gtf_gz=$1 # Annotation in gzipped gtf format
quants_a_tsv=$2      # A quantification file from STAR alignment in tsv format
quants_b_tsv=$3      # B quantification file from STAR alignment in tsv format
out_root="$4_mad"    # Root name for output (e.g. "out" will create "out_mad_qc.txt" and "out_mad_plot.png")

echo "-- Uncompressing annotation..."
annotation_gtf=${annotation_gtf_gz%.gz}
set -x
gunzip $annotation_gtf_gz
set +x

echo "-- Extracting gene ids from annotation..."
set -x
gawk -f extract_gene_ids.awk $annotation_gtf out=srna_gene_ids.txt
set +x
    
echo "-- Generating expression values for '$quants_a_tsv'..."
set -x
gawk -f sum_srna_expression.awk srna_gene_ids.txt $quants_a_tsv out=expr_a.tsv
set +x

echo "-- Generating expression values for '$quants_b_tsv'..."
set -x
gawk -f sum_srna_expression.awk srna_gene_ids.txt $quants_b_tsv out=expr_b.tsv
set +x

echo "-- Runnning MAD.R..."
set -x
Rscript MAD.R expr_a.tsv expr_b.tsv > ${out_root}_qc.txt
mv MAplot.png ${out_root}_plot.png
set +x
echo "-- json? ..."
cat ${out_root}_qc.txt
    
echo "-- The results..."
ls -l ${out_root}*

