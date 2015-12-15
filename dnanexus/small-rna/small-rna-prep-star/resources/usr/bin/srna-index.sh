#!/bin/bash -e

if [ $# -lt 4 ] || [ $# -gt 5 ]; then
    echo "usage v1: srna-index.sh <ref_fasta_gz> <annotation_gtf_gz> <annotation_version> <genome> [<gender>]"
    echo "Converts BAMs from alignments from stranded libraries to bigwig format. Is independent of DX and encodeD."
    exit -1; 
fi
ref_fasta_gz=$1      # Reference genome assembly in gzipped fasta format.
annotation_gtf_gz=$2 # Gene annotation in gzipped gtf format
anno=$3              # Annotation (e.g. 'v24')
genome=$4            # Genome (e.g. 'GRCh38')
if [ $# -eq 5 ]; then
    gender=$5        # Gender. Values: 'female', 'male', 'XX', 'XY' will be included in names.  Otherwise, gender neutral.  
fi

archive_root="${genome}_${anno}_sRNA_starIndex"
if [ "$gender" == "famale" ] || [ "$gender" == "male" ] || [ "$gender" == "XX" ] || [ "$gender" == "XY" ]; then
    archive_root="${genome}_${gender}_${anno}_sRNA_starIndex"
fi
gene_id_root="${genome}_${anno}_gene_ids"
echo "-- Results will be: '${archive_root}.tgz' and '${gene_id_root}.txt'."

echo "-- Unzipping reference files..."
ref_fasta=${ref_fasta_gz%.gz}
ref_root=${ref_fasta%.fasta}
ref_root=${ref_root%.fa}
annotation_gtf=${annotation_gtf_gz%.gz}
gunzip $ref_fasta_gz
gunzip $annotation_gtf_gz

echo "-- Build index for '${genome} and annotation ${anno}'..."
set -x
mkdir out
STAR --runMode genomeGenerate --genomeFastaFiles $ref_fasta --sjdbGTFfile $annotation_gtf \
        --sjdbOverhang 1 --runThreadN 8 --genomeDir out/ --outFileNamePrefix out
set +x
    
echo "-- Build gene ID file from the annotation for later use..."
# Create a geneID file and put it in the index archive for later use:
set -x
gawk -f extract_gene_ids.awk $annotation_gtf out=${gene_id_root}.txt
set +x
#awk '$3=="gene" || substr($14,2,length($14)-3)=="tRNAscan" {g=substr($14,2,length($14)-3); if (g=="miRNA" || g=="snoRNA" || g=="snRNA" || g=="tRNAscan") {print substr($10,2,length($10)-3)} }' \
#    $annotation_gtf | sort > out/smallRNA.geneID
#cp out/smallRNA.geneID ${gene_id_root}.txt

# Attempt to make bamCommentLines.txt, which should be reviewed. NOTE tabs handled by assignment.
echo "-- Create bam header..."
set -x
refComment="@CO\tREFID:$(basename ${ref_root})"
echo -e ${refComment} > out/star_bamCommentLines.txt
echo `cat "out/star_bamCommentLines.txt"`
set +x

echo "-- Create archive file..."
echo "out"
set -x
echo `ls out/`
tar -czf ${archive_root}.tgz out/
set +x

echo "-- The results..."
ls -l ${archive_root}.tgz

