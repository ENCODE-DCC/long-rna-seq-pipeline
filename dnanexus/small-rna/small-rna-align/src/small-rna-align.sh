#!/bin/bash
# small-rna-align 2.0.0

main() {
    # Now in resources/usr/bin
    #echo "* Download and install STAR..."
    #git clone https://github.com/alexdobin/STAR
    #(cd STAR; git checkout tags/STAR_2.4.0k)
    #(cd STAR; make)
    #wget https://github.com/samtools/samtools/archive/0.1.19.tar.gz

    echo "*****"
    echo "* Running: small-rna-align.sh [v2.0.0]"
    echo "* STAR version: "`STAR --version | awk '{print $1}' | cut -d _ -f 2-`
    #echo "* samtools version: "`samtools 2>&1 | grep Version | awk '{print $2}'`
    echo "*****"

    echo "* Value of reads: '$reads'"
    echo "* Value of star_index: '$star_index'"
    echo "* Value of library_id: '$library_id'"
    echo "* Value of nthreads: '$nthreads'"

    echo "* Download files..."
    reads_fn=`dx describe "$reads" --name`
    reads_fn=${reads_fn%.fastq.gz}
    reads_fn=${reads_fn%.fq.gz}
    dx download "$reads" -o "$reads_fn".fastq.gz
    echo "* Read file: '${reads_fn}.fastq.gz'"

    dx download "$star_index" -o star_index.tgz
    tar zxvf star_index.tgz
    # unzips into "out/"

    # Fill in your application code here.

    echo "* Set up headers..."
    libraryComment="@CO\tLIBID:${library_id}"
    echo -e ${libraryComment} > COfile.txt
    cat out/*_bamCommentLines.txt >> COfile.txt
    echo `cat COfile.txt`

    echo "* Map reads..."
    STAR --genomeDir out --readFilesIn ${reads_fn}.fastq.gz --readFilesCommand zcat             \
        --runThreadN ${nthreads} --outFilterMultimapNmax 20 --alignIntronMax 1                  \
        --clip3pAdapterSeq TGGAATTCTC --clip3pAdapterMMp 0.1 --outFilterMismatchNoverLmax 0.05  \
        --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 16  \
        --outSAMheaderCommentFile COfile.txt --outSAMheaderHD @HD VN:1.4 SO:coordinate          \
        --genomeLoad NoSharedMemory --outSAMunmapped Within --outSAMtype BAM SortedByCoordinate \
        --limitBAMsortRAM 60000000000
        
    #echo "* Index genome bam..."
    # Note: no longer making unused index
    mv Aligned.sortedByCoord.out.bam ${reads_fn}_star_genome.bam
    #samtools index ${reads_fn}_star_genome.bam

    # Gather metrics
    meta=`echo \"STAR_log_final\": { `
    #                      Number of input reads |       85357051
    var=`grep "Number of input reads" Log.final.out | awk '{print $6}'`
    var=`echo \"Number of input reads\": $var`
    meta=`echo $meta $var`
    #                  Average input read length |       202
    var=`grep "Average input read length" Log.final.out | awk '{print $6}'`
    var=`echo \"Average input read length\": $var`
    meta=`echo $meta, $var`
    #               Uniquely mapped reads number |       68649148
    var=`grep "Uniquely mapped reads number" Log.final.out | awk '{print $6}'`
    var=`echo \"Uniquely mapped reads number\": $var`
    meta=`echo $meta, $var`
    #                    Uniquely mapped reads % |       80.43%
    var=`grep "Uniquely mapped reads \%" Log.final.out | awk '{print $6}'`
    var=`echo \"Uniquely mapped reads percent\": \"$var\"`
    meta=`echo $meta, $var`
    #                      Average mapped length |       197.40
    var=`grep "Average mapped length" Log.final.out | awk '{print $5}'`
    var=`echo \"Average mapped length\": $var`
    meta=`echo $meta, $var`
    #                   Number of splices: Total |       30116609
    var=`grep "Number of splices\: Total" Log.final.out | awk '{print $6}'`
    var=`echo \"Number of splices total\": $var`
    meta=`echo $meta, $var`
    #        Number of splices: Annotated (sjdb) |       29770156
    var=`grep "Number of splices\: Annotated" Log.final.out | awk '{print $7}'`
    var=`echo \"Number of splices annotated\": $var`
    meta=`echo $meta, $var`
    #                   Number of splices: GT/AG |       29847443
    var=`grep "Number of splices\: GT\/AG" Log.final.out | awk '{print $6}'`
    var=`echo \"Number of splices GT/AG\": $var`
    meta=`echo $meta, $var`
    #                   Number of splices: GC/AG |       207291
    var=`grep "Number of splices\: GC\/AG" Log.final.out | awk '{print $6}'`
    var=`echo \"Number of splices GC/AG\": $var`
    meta=`echo $meta, $var`
    #                   Number of splices: AT/AC |       35287
    var=`grep "Number of splices\: AT\/AC" Log.final.out | awk '{print $6}'`
    var=`echo \"Number of splices AT/AC\": $var`
    meta=`echo $meta, $var`
    #           Number of splices: Non-canonical |       26588
    var=`grep "Number of splices\: Non\-canonical" Log.final.out | awk '{print $6}'`
    var=`echo \"Number of splices non-canonical\": $var`
    meta=`echo $meta, $var`
    #                  Mismatch rate per base, % |       0.62%
    var=`grep "Mismatch rate per base" Log.final.out | awk '{print $7}'`
    var=`echo \"Mismatch rate per base\": \"$var\"`
    meta=`echo $meta, $var`
    #                     Deletion rate per base |       0.01%
    var=`grep "Deletion rate per base" Log.final.out | awk '{print $6}'`
    var=`echo \"Deletion rate per base\": \"$var\"`
    meta=`echo $meta, $var`
    #                    Deletion average length |       1.89
    var=`grep "Deletion average length" Log.final.out | awk '{print $5}'`
    var=`echo \"Deletion average length\": $var`
    meta=`echo $meta, $var`
    #                    Insertion rate per base |       0.02%
    var=`grep "Insertion rate per base" Log.final.out | awk '{print $6}'`
    var=`echo \"Insertion rate per base\": \"$var\"`
    meta=`echo $meta, $var`
    #                   Insertion average length |       1.60
    var=`grep "Insertion average length" Log.final.out | awk '{print $5}'`
    var=`echo \"Insertion average length\": $var`
    meta=`echo $meta, $var`
    #    Number of reads mapped to multiple loci |       3936606
    var=`grep "Number of reads mapped to multiple loci" Log.final.out | awk '{print $9}'`
    var=`echo \"Number of reads mapped to multiple loci\": $var`
    meta=`echo $meta, $var`
    #         % of reads mapped to multiple loci |       4.61%
    var=`grep "\% of reads mapped to multiple loci" Log.final.out | awk '{print $9}'`
    var=`echo \"Percent of reads mapped to multiple loci\": \"$var\"`
    meta=`echo $meta, $var`
    #    Number of reads mapped to too many loci |       10508
    var=`grep "Number of reads mapped to too many loci" Log.final.out | awk '{print $10}'`
    var=`echo \"Number of reads mapped to too many loci\": $var`
    meta=`echo $meta, $var`
    #         % of reads mapped to too many loci |       0.01%
    var=`grep "\% of reads mapped to too many loci" Log.final.out | awk '{print $10}'`
    var=`echo \"Percent of reads mapped to too many loci\": \"$var\"`
    meta=`echo $meta, $var`
    #   % of reads unmapped: too many mismatches |       0.00%
    var=`grep "\% of reads unmapped\: too many mismatches" Log.final.out | awk '{print $9}'`
    var=`echo \"Percent of reads unmapped - too many mismatches\": \"$var\"`
    meta=`echo $meta, $var`
    #             % of reads unmapped: too short |       14.86%
    var=`grep "\% of reads unmapped\: too short" Log.final.out | awk '{print $8}'`
    var=`echo \"Percent reads unmapped - too short\": \"$var\"`
    meta=`echo $meta, $var`
    #                 % of reads unmapped: other |       0.09%
    var=`grep "\% of reads unmapped\: other" Log.final.out | awk '{print $7}'`
    var=`echo \"Percent of reads unmapped - other\": \"$var\"`
    meta=`echo $meta, $var }`

    #mv Aligned.toTranscriptome.out.bam ${reads1_fn}-${reads2_fn}_star_anno.bam
    mv Log.final.out ${reads_fn}_star_Log.final.out

    echo "* Upload results..."
    # NOTE: adding meta 'details' ensures json is valid.  But details are not updatable so rely on QC property
    details=`echo { $meta }`
    genome_bam=$(dx upload ${reads_fn}_star_genome.bam --details "$details" --property QC="$meta" --brief)
    star_log=$(dx upload ${reads_fn}_star_Log.final.out --brief)

    dx-jobutil-add-output genome_bam "$genome_bam" --class=file
    dx-jobutil-add-output star_log "$star_log" --class=file
    dx-jobutil-add-output metadata "$meta" --class=string

    echo "* Finished."
}
