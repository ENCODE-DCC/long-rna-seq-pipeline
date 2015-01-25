# prepare genomes for STAR and RSEM
# input parameters:
STARgenomeDir=$1 #STAR genome directory
RSEMgenomeDir=$2 #RSEM genome directory
fastaGenome=$3   #fasta file(s),  e.g. "male.hg19.fa"
fastaSpikeins=$4 #fasta file with spike-ins, e.g. "spikes.fixed.fasta"
gtf=$5           #all-inclusive gtf file "gencode.v19.annotation.gtf"

# example
# ./STAR_RSEM_prep.sh  /path/to/STARgenome  /path/to/RSEMgenome  male.hg19.fa  spikes.fixed.fasta gencode.v19.annotation_tRNA_spikeins.gtf


## PATHS:
export LD_LIBRARY_PATH=/sonas-hs/gingeras/nlsas_norepl/user/dobin/Software/ZLIB/zlib-1.2.8_installed/lib/
export PATH=$PATH:/sonas-hs/gingeras/nlsas_norepl/user/dobin/Software/RSEM/RSEM-1.2.15/:/sonas-hs/gingeras/nlsas_norepl/user/dobin/STAR/SandBox/STAR/STAR_2.3.1z9/


# RSEM genome
mkdir $RSEMgenomeDir 

### the command below is for RSEM >=1.2.19
### note, that for RSEM < 1.2.19, --no-polyA should be added

RSEMcommand="rsem-prepare-reference --gtf $gtf $fastaGenome","$fastaSpikeins $RSEMgenomeDir/RSEMref"
echo $RSEMcommand
$RSEMcommand


# STAR genome
mkdir $STARgenomeDir
STARcommand="STAR --runThreadN 12 --runMode genomeGenerate --genomeDir $STARgenomeDir --genomeFastaFiles $fastaGenome $fastaSpikeins --sjdbGTFfile $gtf --sjdbOverhang 100 --outFileNamePrefix $STARgenomeDir"
echo $STARcommand
$STARcommand


