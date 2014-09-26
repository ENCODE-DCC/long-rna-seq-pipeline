
/usr/bin/bowtie2-build -f  data/ref/female.hg19.fa,data/refENCFF001RTP.fasta data/female.spikeins.hg19

#- /archive/TGCapps_data/bowtie2_indexes/male.spikeins.hg19

#bowtie2-build -f  /path/to/male.hg19.fa,/path/to/human_spike_ins.fa /archive/TGCapps_data/bowtie2_indexes/male.spikeins.hg19



#Generate transcriptome index file:
#1.       Download genecode.v19.annotation.gtf.gz from http://www.gencodegenes.org/releases/19.html
#2.       Unzip genecode.v19.annotation.gtf.gz
#gunzip gencode.v19.annotation.gtf.gz
#3.       Run tophat with –G option

#- /archive/TGCore_User_Data/ENCODE/ENCODE_ANNOTATIONS/gencode.v19.female.spikeins

/usr/bin/tophat -p 24 -z0 -o output/test-female -a 8 -m 0 --min-intron-length 20 \
  --max-intron-length 1000000 --read-edit-dist 4 --read-mismatches 4 -g 20 \
  --no-discordant --no-mixed --library-type fr-firststrand -G data/gencode.v19.annotation.gtf \
  --transcriptome-index data/gencode.v19.female.spikeins data/female.spikeins.hg19 data/Read1.fastq.gz data/Read2.fastq.gz

#- /archive/TGCore_User_Data/ENCODE/ENCODE_ANNOTATIONS/gencode.v19.male.spikeins

#/usr/bin/tophat -p 24 -z0 -o /output/folder/tophat/test2-male -a 8 -m 0 --min-intron-length 20 --max-intron-length 1000000 --read-edit-dist 4 --read-mismatches 4 -g 20 --no-discordant --no-mixed --library-type fr-firststrand -G /archive/TGCore_User_Data/ENCODE/ENCODE_ANNOTATIONS/gencode.v19.annotation.gtf --transcriptome-index /archive/TGCore_User_Data/ENCODE/ENCODE_ANNOTATIONS/gencode.v19.male.spikeins /archive/TGCapps_data/bowtie2_indexes/male.spikeins.hg19 /path/to/Read1.fastq.gz /path/to/Read2.fastq.gz



#Generate a chromosome sizes file (should use male’s bam file):
#- /archive/TGCore_User_Data/ENCODE/ENCODE_ANNOTATIONS/chromosizes_spikeIns.txt

/usr/binsamtools view -H /output/folder/tophat/test2-male/accepted_hits.bam | awk '/@SQ/ {OFS="\t"; gsub("SN:", "", $2); gsub("LN:", "", $3); print $2, $3}' > /archive/TGCore_User_Data/ENCODE/ENCODE_ANNOTATIONS/chromosizes_spikeIns.txt

