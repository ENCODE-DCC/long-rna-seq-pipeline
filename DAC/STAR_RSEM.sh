#!/bin/bash

# STAR mapping / RSEM quantification pipeline
# usage: from an empty working directory, run
# ./STAR_RSEM.sh (read1) (read2 or "") (STARgenomeDir) (RSEMrefDir) (dataType) (nThreadsSTAR) (nThreadsRSEM)

# input: gzipped fastq file read1 [read2 for paired-end] 
#        STAR genome directory, RSEM reference directory - prepared with STAR_RSEM_prep.sh script
read1=$1 #gzipped fastq file for read1
read2=$2 #gzipped fastq file for read1, use "" if single-end
STARgenomeDir=$3 
RSEMrefDir=$4
dataType=$5 # RNA-seq type, possible values: str_SE str_PE unstr_SE unstr_PE
nThreadsSTAR=$6 # number of threads for STAR
nThreadsRSEM=$7 # number of threads for RSEM

# output: all in the working directory, fixed names
# Aligned.sortedByCoord.out.bam                 # alignments, standard sorted BAM, agreed upon formatting
# Log.final.out                                 # mapping statistics to be used for QC, text, STAR formatting
# Quant.genes.results                           # RSEM gene quantifications, tab separated text, RSEM formatting
# Quant.isoforms.results                        # RSEM transcript quantifications, tab separated text, RSEM formatting
# Quant.pdf                                     # RSEM diagnostic plots
# Signal.{Unique,UniqueMultiple}.strand{+,-}.bw # 4 bigWig files for stranded data
# Signal.{Unique,UniqueMultiple}.unstranded.bw  # 2 bigWig files for unstranded data

# executables
STAR=STAR                             
RSEM=rsem-calculate-expression        
bedGraphToBigWig=bedGraphToBigWig              


# STAR parameters: common
STARparCommon=" --genomeDir $STARgenomeDir  --readFilesIn $read1 $read2   --outSAMunmapped Within --outFilterType BySJout \
 --outSAMattributes NH HI AS NM MD    --outFilterMultimapNmax 20   --outFilterMismatchNmax 999   \
 --outFilterMismatchNoverReadLmax 0.04   --alignIntronMin 20   --alignIntronMax 1000000   --alignMatesGapMax 1000000   \
 --alignSJoverhangMin 8   --alignSJDBoverhangMin 1 --sjdbScore 1 --readFilesCommand zcat"

# STAR parameters: run-time, controlled by DCC
STARparRun=" --runThreadN $nThreadsSTAR --genomeLoad LoadAndKeep  --limitBAMsortRAM 10000000000"

# STAR parameters: type of BAM output: quantification or sorted BAM or both
#     OPTION: sorted BAM output
## STARparBAM="--outSAMtype BAM SortedByCoordinate"
#     OPTION: transcritomic BAM for quantification
## STARparBAM="--outSAMtype None --quantMode TranscriptomeSAM"
#     OPTION: both
STARparBAM="--outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM"


# STAR parameters: strandedness, affects bedGraph (wiggle) files and XS tag in BAM 

case "$dataType" in
str_SE|str_PE)
      #OPTION: stranded data
      STARparStrand=""
      STARparWig="--outWigStrand Stranded"
      ;;
      #OPTION: unstranded data
unstr_SE|unstr_PE)
      STARparStrand="--outSAMstrandField intronMotif"
      STARparWig="--outWigStrand Unstranded"
      ;;
esac

# STAR parameters: metadata
STARparsMeta="--outSAMheaderCommentFile commentsENCODElong.txt --outSAMheaderHD @HD VN:1.4 SO:coordinate"

## not needed ## --outSAMheaderPG @PG ID:Samtools PN:Samtools CL:"$samtoolsCommand" PP:STAR VN:0.1.18"

# ENCODE metadata BAM comments
echo -e '@CO\tLIBID:ENCLB175ZZZ
@CO\tREFID:ENCFF001RGS
@CO\tANNID:gencode.v19.annotation.gtf.gz
@CO\tSPIKEID:ENCFF001RTP VN:Ambion-ERCC Mix, Cat no. 445670' > commentsENCODElong.txt

###### STAR command
echo $STAR $STARparCommon $STARparRun $STARparBAM $STARparStrand $STARparsMeta
$STAR $STARparCommon $STARparRun $STARparBAM $STARparStrand $STARparsMeta

###### bedGraph generation, now decoupled from STAR alignment step
# working subdirectory for this STAR run
mkdir Signal

echo $STAR --runMode inputAlignmentsFromBAM   --inputBAMfile Aligned.sortedByCoord.out.bam --outWigType bedGraph $STARparWig --outFileNamePrefix ./Signal/ --outWigReferencesPrefix chr
$STAR --runMode inputAlignmentsFromBAM   --inputBAMfile Aligned.sortedByCoord.out.bam --outWigType bedGraph $STARparWig --outFileNamePrefix ./Signal/ --outWigReferencesPrefix chr

# move the signal files from the subdirectory
mv Signal/Signal*bg .




###### bigWig conversion commands
# exclude spikeins
grep ^chr $STARgenomeDir/chrNameLength.txt > chrNL.txt

case "$dataType" in
str_SE|str_PE)
      # stranded data
      str[1]=-; str[2]=+;
      for istr in 1 2
      do
      for imult in Unique UniqueMultiple
      do
          grep ^chr Signal.$imult.str$istr.out.bg > sig.tmp
          $bedGraphToBigWig sig.tmp  chrNL.txt Signal.$imult.strand${str[istr]}.bw
      done
      done
      ;;
unstr_SE|unstr_PE)
      # unstranded data
      for imult in Unique UniqueMultiple
      do
          grep ^chr Signal.$imult.str1.out.bg > sig.tmp
          $bedGraphToBigWig sig.tmp chrNL.txt  Signal.$imult.unstranded.bw
      done
      ;;
esac




######### RSEM

#### prepare for RSEM: sort transcriptome BAM to ensure the order of the reads, to make RSEM output (not pme) deterministic
trBAMsortRAM=60G

mv Aligned.toTranscriptome.out.bam Tr.bam 

case "$dataType" in
str_SE|unstr_SE)
      # single-end data
      cat <( samtools view -H Tr.bam ) <( samtools view -@ $nThreadsRSEM Tr.bam | sort -S $trBAMsortRAM -T ./ ) | samtools view -@ $nThreadsRSEM -bS - > Aligned.toTranscriptome.out.bam
      ;;
str_PE|unstr_PE)
      # paired-end data, merge mates into one line before sorting, and un-merge after sorting
      cat <( samtools view -H Tr.bam ) <( samtools view -@ $nThreadsRSEM Tr.bam | awk '{printf "%s", $0 " "; getline; print}' | sort -S $trBAMsortRAM -T ./ | tr ' ' '\n' ) | samtools view -@ $nThreadsRSEM -bS - > Aligned.toTranscriptome.out.bam
      ;;
esac

'rm' Tr.bam


# RSEM parameters: common
RSEMparCommon="--bam --estimate-rspd  --calc-ci --no-bam-output --seed 12345"

# RSEM parameters: run-time, number of threads and RAM in MB
RSEMparRun=" -p $nThreadsRSEM --ci-memory 30000 "

# RSEM parameters: data type dependent

case "$dataType" in
str_SE)
      #OPTION: stranded single end
      RSEMparType="--forward-prob 0"
      ;;
str_PE)
      #OPTION: stranded paired end
      RSEMparType="--paired-end --forward-prob 0"
      ;;
unstr_SE)
      #OPTION: unstranded single end
      RSEMparType=""
      ;;
unstr_PE)
      #OPTION: unstranded paired end
      RSEMparType="--paired-end"
      ;;
esac


###### RSEM command
echo $RSEM $RSEMparCommon $RSEMparRun $RSEMparType Aligned.toTranscriptome.out.bam $RSEMrefDir Quant >& Log.rsem
$RSEM $RSEMparCommon $RSEMparRun $RSEMparType Aligned.toTranscriptome.out.bam $RSEMrefDir Quant >& Log.rsem

###### RSEM diagnostic plot creation
# Notes:
# 1. rsem-plot-model requires R (and the Rscript executable)
# 2. This command produces the file Quant.pdf, which contains multiple plots
echo rsem-plot-model Quant Quant.pdf
rsem-plot-model Quant Quant.pdf
