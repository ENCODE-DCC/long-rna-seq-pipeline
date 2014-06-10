# prepare genomes for STAR and RSEM
input parameters:
# $1=fasta file(s),  e.g. "male.hg19.fa"
# $2=fasta file with spike-ins, e.g. "spikes.fixed.fasta"
# $3=gtf file, e.g. "gencode.v19.annotation.gtf"
# $4=STAR genome directory
# $5=RSEM genome director

# example
./STAR_RSEM_prep.sh male.hg19.fa  spikes.fixed.fasta gencode.v19.annotation.gtf /path/to/STARgenome /path/to/RSEMgenome


mkdir $4
cd $4
STAR --runMode genomeGenerate --runThreadN 12 --genomeDir ./ --genomeFastaFiles  $1 $2 --sjdbGTFfile $3 --sjdbOverhang 100

mkdir $5
cd $5
rsem-prepare-reference --no-polyA --no-bowtie --no-ntog --gtf $3 $1 RSEMref
