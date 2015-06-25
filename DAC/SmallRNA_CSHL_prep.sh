# prepare genomes for STAR and RSEM
# run from the directory intended to be STAR genome directory
# input parameters:
fastaGenome=$1   #fasta file(s),  e.g. "male.hg19.fa"
gtf=$2           #all-inclusive gtf file "gencode.v19.annotation.gtf"

# example
# ./STAR_RSEM_prep.sh male.hg19.fa  gencode.v19.annotation_tRNA_spikeins.gtf


## PATHS:
export LD_LIBRARY_PATH=/sonas-hs/gingeras/nlsas_norepl/user/dobin/Software/ZLIB/zlib-1.2.8_installed/lib/
export PATH=$PATH:/sonas-hs/gingeras/nlsas_norepl/user/dobin/Software/RSEM/RSEM-1.2.15/:/sonas-hs/gingeras/nlsas_norepl/user/dobin/STAR/SandBox/STAR/STAR_2.3.1z9/

# small RNAs gene ID for quantifications and QC
awk '$3=="gene" || substr($14,2,length($14)-3)=="tRNAscan" {g=substr($14,2,length($14)-3); if (g=="miRNA" || g=="snoRNA" || g=="snRNA" || g=="tRNAscan") {print substr($10,2,length($10)-3)} } ' $gtf > smallRNA.geneID

# STAR genome
STARcommand="STAR --runThreadN 12 --runMode genomeGenerate --genomeDir ./ --genomeFastaFiles $fastaGenome --sjdbGTFfile $gtf --sjdbOverhang 1 "
echo $STARcommand
$STARcommand



