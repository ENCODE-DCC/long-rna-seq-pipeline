# STAR mapping / RSEM quantification pipeline
# usage: from an empty working directory, run
# ./STAR_RSEM.sh <read1> <read2 or ""> <STARgenomeDir> <RSEMrefDir>

# input: gzipped fastq file read1 [read2 for paired-end] 
#        STAR genome directory, RSEM reference directory - prepared with STAR_RSEM_prep.sh script
read1=$1 #gzipped fastq file for read1
read2=$2 #gzipped fastq file for read1, use "" if single-end
STARgenomeDir=$3 
RSEMrefDir=$4

# output: all in the working directory, fixed names
# Aligned.sortedByCoord.out.bam                 # alignments, standard sorted BAM, agreed upon formatting
# Log.final.out                                 # mapping statistics to be used for QC, text, STAR formatting
# Quant.genes.results                           # RSEM gene quantifications, tab separated text, RSEM formatting
# Quant.isoforms.results                        # RSEM transcript quantifications, tab separated text, RSEM formatting
# Signal.{Unique,UniqueMultiple}.strand{+,-}.bw # 4 bigWig files for stranded data
# Signal.{Unique,UniqueMultiple}.unstranded.bw  # 4 bigWig files for stranded data

# executables
STAR=./STAR                             # version: 2.4.0,   GitHub link:
RSEM=./rsem-1.2.7/rsem-calculate-expression        # version: x.x.x,   GitHub link:
wigToBigWig=./wigToBigWig               # version: v4,      GitHub link:


# STAR parameters: common
STARparCommon="--genomeDir $STARgenomeDir  --readFilesIn $read1 $read2 \
 --outSAMattributes NH HI AS NM MD    --outFilterMultimapNmax 20   --outFilterMismatchNmax 999   \
 --outFilterMismatchNoverLmax 0.04   --alignIntronMin 20   --alignIntronMax 1000000   --alignMatesGapMax 1000000   \
 --alignSJoverhangMin 8   --alignSJDBoverhangMin 1 --readFilesCommand zcat --outWigType bedGraph"

# STAR parameters: run-time, controlled by DCC
STARparRun="--runThreadN 12 --genomeLoad LoadAndKeep"

# STAR parameters: type of BAM output: quantification or sorted BAM or both
#     OPTION: sorted BAM output
## STARparBAM="--outSAMtype BAM SortedByCoordinate"
#     OPTION: transcritomic BAM for quantification
## STARparBAM="--outSAMtype None --quantMode TranscriptomeSAM"
#     OPTION: both
STARparBAM="--outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM"

# STAR parameters: strandedness, affects bedGraph (wiggle) files and XS tag in BAM 
#     OPTION: stranded data
# STARparStrand="--outWigStrand Stranded"
#     OPTION: unstranded data
STARparStrand="--outWigStrand Unstranded --outSAMstrandField intronMotif"

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



# RSEM parameters: common
RSEMparCommon="--bam --estimate-rspd  --calc-ci --no-bam-output"

# RSEM parameters: run-time, number of threads and RAM in MB
RSEMparRun="-p 12 --ci-memory 30000"

# RSEM parameters: data type dependent
#    OPTION: stranded paired end
# RSEMparType="--paired-end --forward-prob 0"
#    OPTION: unstranded single end
# RSEMparType=""
#    OPTION: unstranded paired end
RSEMparType="--paired-end"
#    OPTION: stranded single end
# RSEMparType="--forward-prob 0"



###### RSEM command
echo $RSEM $RSEMparCommon $RSEMparRun $RSEMparType Aligned.toTranscriptome.out.bam $RSEMrefDir Quant
$RSEM $RSEMparCommon $RSEMparRun $RSEMparType Aligned.toTranscriptome.out.bam $RSEMrefDir Quant >& Log.rsem



# ??? we may need a simple script here to filter RSEM results for unwanted genes

###### bigWig conversion commands: stranded data
#str[1]=-; str[2]=+;
#for istr in 1 2
#do
#for imult in Unique UniqueMultiple
#do
#    $wigToBigWig Signal.$imult.str$istr.out.bg $STARgenomeDir/chrNameLength.txt Signal.$imult.strand${str[istr]}.bw
#done
#done
###### bigWig conversion commands: unstranded data
for imult in Unique UniqueMultiple
do
    $wigToBigWig Signal.$imult.str1.out.bg $STARgenomeDir/chrNameLength.txt Signal.$imult.unstranded.bw
done
