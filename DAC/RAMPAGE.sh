#!/bin/bash

# STAR mapping / RSEM quantification pipeline
# usage: from an empty working directory, run
# ./STAR_RSEM.sh (read1) (read2 or "") (STARgenomeDir) (RSEMrefDir) (dataType) (nThreadsSTAR)

# input: gzipped fastq file read1 [read2 for paired-end] 
#        STAR genome directory, RSEM reference directory - prepared with STAR_RSEM_prep.sh script
read1=$1 #gzipped fastq file for read1
read2=$2 #gzipped fastq file for read1, use "" if single-end
STARgenomeDir=$3 
nThreadsSTAR=$4 # number of threads for STAR

# output: all in the working directory, fixed names
# markdup.Processed.out.bam                     # alignments, standard sorted BAM, agreed upon formatting, PCR duplicates marked
# Log.final.out                                 # mapping statistics to be used for QC, text, STAR formatting
# Signal.{Unique,UniqueMultiple}.strand{+,-}.bw # 4 bigWig files for stranded data

# executables
STAR=STAR                             
bedGraphToBigWig=bedGraphToBigWig              

# STAR parameters: common
STARparCommon="--genomeDir $STARgenomeDir  --readFilesIn $read1 $read2   --outSAMunmapped Within --outFilterType BySJout \
 --outSAMattributes NH HI AS NM MD    --outFilterMultimapNmax 500   --outFilterMismatchNmax 999   \
 --outFilterMismatchNoverReadLmax 0.04   --alignIntronMin 20   --alignIntronMax 1000000   --alignMatesGapMax 1000000   \
 --alignSJoverhangMin 8   --alignSJDBoverhangMin 1 --readFilesCommand zcat  --outSAMtype BAM SortedByCoordinate  \
 --outFilterScoreMinOverLread 0.85  --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
 --clip5pNbases 6   15  --seedSearchStartLmax 30 "

# STAR parameters: run-time, controlled by DCC
STARsortRAM=" --limitBAMsortRAM 30000000000 "
STARparRun=" --runThreadN $nThreadsSTAR --genomeLoad LoadAndKeep $STARsortRAM "


# STAR parameters: metadata
STARparsMeta="--outSAMheaderCommentFile commentsENCODElong.txt --outSAMheaderHD @HD VN:1.4 SO:coordinate"

# ENCODE metadata BAM comments
echo -e '@CO\tLIBID:ENCLB175ZZZ
@CO\tREFID:ENCFF001RGS
@CO\tANNID:gencode.v19.annotation.gtf.gz
@CO\tSPIKEID:ENCFF001RTP VN:Ambion-ERCC Mix, Cat no. 445670' > commentsENCODElong.txt

###### STAR command
echo $STAR $STARparCommon $STARparRun $STARparsMeta
     $STAR $STARparCommon $STARparRun $STARparsMeta

#####################################################################################
### mark PCR duplicates:
#####################################################################################

echo $STAR --inputBAMfile Aligned.sortedByCoord.out.bam --bamRemoveDuplicatesType UniqueIdentical --runMode inputAlignmentsFromBAM \
--bamRemoveDuplicatesMate2basesN 15 --outFileNamePrefix markdup.  $STARsortRAM
     $STAR --inputBAMfile Aligned.sortedByCoord.out.bam --bamRemoveDuplicatesType UniqueIdentical --runMode inputAlignmentsFromBAM \
--bamRemoveDuplicatesMate2basesN 15 --outFileNamePrefix markdup.  $STARsortRAM


#####################################################################################
### signal tracks
#####################################################################################

mkdir Signal
cd Signal

## Generate 5' end signal file:
echo $STAR --inputBAMfile ../markdup.Processed.out.bam --runMode inputAlignmentsFromBAM  --outWigType bedGraph read1_5p --outFileNamePrefix read1_5p. --outWigReferencesPrefix chr
     $STAR --inputBAMfile ../markdup.Processed.out.bam --runMode inputAlignmentsFromBAM  --outWigType bedGraph read1_5p --outFileNamePrefix read1_5p. --outWigReferencesPrefix chr

## bigWig conversion commands
# exclude spikeins
grep ^chr $STARgenomeDir/chrNameLength.txt > chrNL.txt

# stranded data only
str[1]=+; str[2]=-;
for istr in 1 2
do
for imult in Unique UniqueMultiple
    do
      $bedGraphToBigWig read1_5p.Signal.$imult.str$istr.out.bg  chrNL.txt read1_5p.$imult.strand${str[istr]}.bw
    done
done

mv *.bw ..

cd ..


#####################################################################################
### call RAMPAGE peaks
#####################################################################################

