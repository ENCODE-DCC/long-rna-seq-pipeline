#!/bin/bash
# align-star-pe.sh

script_name="align-star-pe.sh"
script_ver="2.0.1"

main() {
    # Now in resources/usr/bin
    #echo "* Download and install STAR..."
    #git clone https://github.com/alexdobin/STAR
    #(cd STAR; git checkout tags/STAR_2.4.0k)
    #(cd STAR; make)
    #wget https://github.com/samtools/samtools/archive/0.1.19.tar.gz

    # If available, will print tool versions to stderr and json string to stdout
    versions=''
    if [ -f /usr/bin/tool_versions.py ]; then 
        versions=`tool_versions.py --applet $script_name --appver $script_ver`
    fi
    #echo "*****"
    #echo "* Running: align-star-pe.sh [v2.0.0]"
    #echo "* STAR version: "`STAR --version | awk '{print $1}' | cut -d _ -f 2-`
    #echo "* samtools version: "`samtools 2>&1 | grep Version | awk '{print $2}'`
    #echo "*****"

    echo "* Value of reads_1: '$reads_1'"
    echo "* Value of reads_2: '$reads_2'"
    echo "* Value of star_index: '$star_index'"
    echo "* Value of library_id: '$library_id'"
    echo "* Number of threads (default 8): '$nthreads'"

    echo "* Download files..."
    reads1_fn=`dx describe "$reads_1" --name`
    reads1_fn=${reads1_fn%.fastq.gz}
    reads1_fn=${reads1_fn%.fq.gz}
    dx download "$reads_1" -o "$reads1_fn".fastq.gz
    reads2_fn=`dx describe "$reads_2" --name`
    reads2_fn=${reads2_fn%.fastq.gz}
    reads2_fn=${reads2_fn%.fq.gz}
    dx download "$reads_2" -o "$reads2_fn".fastq.gz
    echo "* Read files: '${reads1_fn}.fastq.gz' '${reads2_fn}.fastq.gz'"

    dx download "$star_index" -o star_index.tgz
    tar zxvf star_index.tgz
    # unzips into "out/"

    # Fill in your application code here.

    echo "* Set up headers..."
    set -x
    libraryComment="@CO\tLIBID:${library_id}"
    echo -e ${libraryComment} > COfile.txt
    cat out/*_bamCommentLines.txt >> COfile.txt
    set +x
    echo `cat COfile.txt`

    echo "* Map reads..."
    set -x
    STAR --genomeDir out --readFilesIn ${reads1_fn}.fastq.gz ${reads2_fn}.fastq.gz  \
        --readFilesCommand zcat --runThreadN ${nthreads} --genomeLoad NoSharedMemory \
        --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1    \
        --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04              \
        --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000         \
        --outSAMheaderCommentFile COfile.txt --outSAMheaderHD @HD VN:1.4 SO:coordinate   \
        --outSAMunmapped Within --outFilterType BySJout --outSAMattributes NH HI AS NM MD \
        --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM --sjdbScore 1     \
        --limitBAMsortRAM 60000000000

    #echo "* Index genome bam..."
    # Note: no longer making unused index
    mv Aligned.sortedByCoord.out.bam ${reads1_fn}-${reads2_fn}_star_genome.bam
    #samtools index ${reads1_fn}-${reads2_fn}_star_genome.bam
    set +x

    echo "* Sorting annotation bam..."
    set -x
    cat <( samtools view -H Aligned.toTranscriptome.out.bam ) \
        <( samtools view -@ ${nthreads} Aligned.toTranscriptome.out.bam | \
            awk '{printf $0 " "; getline; print}' | \
            sort -S 60G -T ./ | tr ' ' '\n' ) | \
        samtools view -@ ${nthreads} -bS - > ${reads1_fn}-${reads2_fn}_star_anno.bam
    set +x
    
    echo "* Prepare metadata..."
    meta=''
    if [ -f /usr/bin/qc_metrics.py ]; then
        meta=`qc_metrics.py -n STAR_log_final -f Log.final.out`
    fi
    ## Gather metrics
    #meta=`echo \"STAR_log_final\": { `
    ##                      Number of input reads |       85357051
    #var=`grep "Number of input reads" Log.final.out | awk '{print $6}'`
    #var=`echo \"Number of input reads\": $var`
    #meta=`echo $meta $var`
    ##                  Average input read length |       202
    #var=`grep "Average input read length" Log.final.out | awk '{print $6}'`
    #var=`echo \"Average input read length\": $var`
    #meta=`echo $meta, $var`
    ##               Uniquely mapped reads number |       68649148
    #var=`grep "Uniquely mapped reads number" Log.final.out | awk '{print $6}'`
    #var=`echo \"Uniquely mapped reads number\": $var`
    #meta=`echo $meta, $var`
    ##                    Uniquely mapped reads % |       80.43%
    #var=`grep "Uniquely mapped reads \%" Log.final.out | awk '{print $6}'`
    #var=`echo \"Uniquely mapped reads percent\": \"$var\"`
    #meta=`echo $meta, $var`
    ##                      Average mapped length |       197.40
    #var=`grep "Average mapped length" Log.final.out | awk '{print $5}'`
    #var=`echo \"Average mapped length\": $var`
    #meta=`echo $meta, $var`
    ##                   Number of splices: Total |       30116609
    #var=`grep "Number of splices\: Total" Log.final.out | awk '{print $6}'`
    #var=`echo \"Number of splices total\": $var`
    #meta=`echo $meta, $var`
    ##        Number of splices: Annotated (sjdb) |       29770156
    #var=`grep "Number of splices\: Annotated" Log.final.out | awk '{print $7}'`
    #var=`echo \"Number of splices annotated\": $var`
    #meta=`echo $meta, $var`
    ##                   Number of splices: GT/AG |       29847443
    #var=`grep "Number of splices\: GT\/AG" Log.final.out | awk '{print $6}'`
    #var=`echo \"Number of splices GT/AG\": $var`
    #meta=`echo $meta, $var`
    ##                   Number of splices: GC/AG |       207291
    #var=`grep "Number of splices\: GC\/AG" Log.final.out | awk '{print $6}'`
    #var=`echo \"Number of splices GC/AG\": $var`
    #meta=`echo $meta, $var`
    ##                   Number of splices: AT/AC |       35287
    #var=`grep "Number of splices\: AT\/AC" Log.final.out | awk '{print $6}'`
    #var=`echo \"Number of splices AT/AC\": $var`
    #meta=`echo $meta, $var`
    ##           Number of splices: Non-canonical |       26588
    #var=`grep "Number of splices\: Non\-canonical" Log.final.out | awk '{print $6}'`
    #var=`echo \"Number of splices non-canonical\": $var`
    #meta=`echo $meta, $var`
    ##                  Mismatch rate per base, % |       0.62%
    #var=`grep "Mismatch rate per base" Log.final.out | awk '{print $7}'`
    #var=`echo \"Mismatch rate per base\": \"$var\"`
    #meta=`echo $meta, $var`
    ##                     Deletion rate per base |       0.01%
    #var=`grep "Deletion rate per base" Log.final.out | awk '{print $6}'`
    #var=`echo \"Deletion rate per base\": \"$var\"`
    #meta=`echo $meta, $var`
    ##                    Deletion average length |       1.89
    #var=`grep "Deletion average length" Log.final.out | awk '{print $5}'`
    #var=`echo \"Deletion average length\": $var`
    #meta=`echo $meta, $var`
    ##                    Insertion rate per base |       0.02%
    #var=`grep "Insertion rate per base" Log.final.out | awk '{print $6}'`
    #var=`echo \"Insertion rate per base\": \"$var\"`
    #meta=`echo $meta, $var`
    ##                   Insertion average length |       1.60
    #var=`grep "Insertion average length" Log.final.out | awk '{print $5}'`
    #var=`echo \"Insertion average length\": $var`
    #meta=`echo $meta, $var`
    ##    Number of reads mapped to multiple loci |       3936606
    #var=`grep "Number of reads mapped to multiple loci" Log.final.out | awk '{print $9}'`
    #var=`echo \"Number of reads mapped to multiple loci\": $var`
    #meta=`echo $meta, $var`
    ##         % of reads mapped to multiple loci |       4.61%
    #var=`grep "\% of reads mapped to multiple loci" Log.final.out | awk '{print $9}'`
    #var=`echo \"Percent of reads mapped to multiple loci\": \"$var\"`
    #meta=`echo $meta, $var`
    ##    Number of reads mapped to too many loci |       10508
    #var=`grep "Number of reads mapped to too many loci" Log.final.out | awk '{print $10}'`
    #var=`echo \"Number of reads mapped to too many loci\": $var`
    #meta=`echo $meta, $var`
    ##         % of reads mapped to too many loci |       0.01%
    #var=`grep "\% of reads mapped to too many loci" Log.final.out | awk '{print $10}'`
    #var=`echo \"Percent of reads mapped to too many loci\": \"$var\"`
    #meta=`echo $meta, $var`
    ##   % of reads unmapped: too many mismatches |       0.00%
    #var=`grep "\% of reads unmapped\: too many mismatches" Log.final.out | awk '{print $9}'`
    #var=`echo \"Percent of reads unmapped - too many mismatches\": \"$var\"`
    #meta=`echo $meta, $var`
    ##             % of reads unmapped: too short |       14.86%
    #var=`grep "\% of reads unmapped\: too short" Log.final.out | awk '{print $8}'`
    #var=`echo \"Percent reads unmapped - too short\": \"$var\"`
    #meta=`echo $meta, $var`
    ##                 % of reads unmapped: other |       0.09%
    #var=`grep "\% of reads unmapped\: other" Log.final.out | awk '{print $7}'`
    #var=`echo \"Percent of reads unmapped - other\": \"$var\"`
    #meta=`echo $meta, $var }`
    #echo "metadata = <" ${meta} ">"

    #mv Aligned.toTranscriptome.out.bam ${reads1_fn}-${reads2_fn}_star_anno.bam
    set -x
    mv Log.final.out ${reads1_fn}-${reads2_fn}_star_Log.final.out
    set +x

    echo "* Upload results..."
    # NOTE: adding meta 'details' ensures json is valid.  But details are not updatable so rely on QC property
    star_genome_bam=$(dx upload ${reads1_fn}-${reads2_fn}_star_genome.bam --details="{ $meta }" --property QC="{ $meta }" --property SW="$versions" --brief)
    star_anno_bam=$(dx upload ${reads1_fn}-${reads2_fn}_star_anno.bam     --details="{ $meta }" --property QC="{ $meta }" --property SW="$versions" --brief)
    star_log=$(dx upload ${reads1_fn}-${reads2_fn}_star_Log.final.out --property SW="$versions" --brief)

    dx-jobutil-add-output star_genome_bam "$star_genome_bam" --class=file
    dx-jobutil-add-output star_anno_bam "$star_anno_bam" --class=file
    dx-jobutil-add-output star_log "$star_log" --class=file
    dx-jobutil-add-output metadata "$meta" --class=string

    echo "* Finished."
}
