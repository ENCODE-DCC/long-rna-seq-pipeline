#!/bin/bash

# CSHL small RNA pipeline: STAR mapping to the genome, bigWig generation
# usage: from an empty working directory, run
# ./STAR_RSEM.sh (read1) (STARgenomeDir) (nThreadsSTAR)

# input: gzipped fastq file read1 [read2 for paired-end] 
#        STAR genome directory, RSEM reference directory - prepared with STAR_RSEM_prep.sh script
read1=$1 #gzipped fastq file for read1
STARgenomeDir=$2
nThreadsSTAR=$3 # number of threads for STAR

# output: all in the working directory, fixed names
# Aligned.sortedByCoord.out.bam                 # alignments, standard sorted BAM, agreed upon formatting, PCR duplicates marked
# Log.final.out                                 # mapping statistics to be used for QC, text, STAR formatting
# Signal.{Unique,UniqueMultiple}.strand{+,-}.bw # 4 bigWig files for stranded data

# executables
STAR=STAR                             
bedGraphToBigWig=bedGraphToBigWig              

# STAR parameters: common
STARparCommon=" --genomeDir $STARgenomeDir  --readFilesIn $read1   --outSAMunmapped Within \
                --readFilesCommand zcat  --outSAMtype BAM SortedByCoordinate  --quantMode GeneCounts \
                --outFilterMultimapNmax 20 --clip3pAdapterSeq TGGAATTCTC --clip3pAdapterMMp 0.1 --outFilterMismatchNoverLmax 0.03 \
                --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 16 --alignSJDBoverhangMin 1000 --alignIntronMax 1 "


# STAR parameters: run-time, controlled by DCC
STARsortRAM=" --limitBAMsortRAM 30000000000 "
STARparRun=" --runThreadN $nThreadsSTAR --genomeLoad LoadAndKeep $STARsortRAM "


# STAR parameters: metadata
STARparsMeta="--outSAMheaderCommentFile commentsENCODElong.txt --outSAMheaderHD @HD VN:1.4 SO:coordinate"

# ENCODE metadata BAM comments
echo -e '@CO\tLIBID:ENCLB175ZZZ
@CO\tREFID:ENCFF001RGS' > commentsENCODElong.txt

###### STAR command
echo $STAR $STARparCommon $STARparRun $STARparsMeta
     $STAR $STARparCommon $STARparRun $STARparsMeta


#####################################################################################
### reads per gene for small RNA
#####################################################################################
# pre-processing command for annotations
awk 'BEGIN {OFS="\t"} {if (ARGIND==1) {G[$1]=1} else {if (substr($1,1,2)!="N_") {t+=$3} else {print $1,$3}; if ($1 in G) {print $1,".",".",".",".",".",$3 >"smallRNA.expr"; s+=$3 }}} END {print "N_notSmallRNA", t-s; print "N_smallRNA",s}' $STARgenomeDir/smallRNA.geneID ReadsPerGene.out.tab > smallRNA.info

#####################################################################################
### signal tracks
#####################################################################################

mkdir Signal
cd Signal

## Generate signal file:
echo $STAR --inputBAMfile ../Aligned.sortedByCoord.out.bam --runMode inputAlignmentsFromBAM  --outWigType bedGraph  --outWigReferencesPrefix chr
     $STAR --inputBAMfile ../Aligned.sortedByCoord.out.bam --runMode inputAlignmentsFromBAM  --outWigType bedGraph  --outWigReferencesPrefix chr

## bigWig conversion commands
# exclude spikeins
grep ^chr $STARgenomeDir/chrNameLength.txt > chrNL.txt

# stranded data only
str[1]=+; str[2]=-;
for istr in 1 2
do
for imult in Unique UniqueMultiple
    do
      $bedGraphToBigWig Signal.$imult.str$istr.out.bg  chrNL.txt Signal.$imult.strand${str[istr]}.bw
    done
done

mv *.bw ..

cd ..

