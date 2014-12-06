# prepare genomes for STAR and RSEM
# input parameters:
STARgenomeDir=$1 #STAR genome directory
fastaGenome=$2   #fasta file(s),  e.g. "male.hg19.fa"

# example
# ./SmallRNA_CSHL_prep.sh  /path/to/STARgenome  male.hg19.fa


## PATHS:
export LD_LIBRARY_PATH=/sonas-hs/gingeras/nlsas_norepl/user/dobin/Software/ZLIB/zlib-1.2.8_installed/lib/
export PATH=$PATH:/sonas-hs/gingeras/nlsas_norepl/user/dobin/Software/RSEM/RSEM-1.2.15/:/sonas-hs/gingeras/nlsas_norepl/user/dobin/STAR/SandBox/STAR/STAR_2.3.1z9/


# STAR genome
mkdir $STARgenomeDir
STARcommand="STAR --runThreadN 12 --runMode genomeGenerate --genomeDir $STARgenomeDir --genomeFastaFiles $fastaGenome --outFileNamePrefix $STARgenomeDir"
echo $STARcommand
$STARcommand


