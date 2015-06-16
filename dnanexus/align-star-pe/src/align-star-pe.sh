#!/bin/bash
# align-star-pe.sh

script_name="align-star-pe.sh"
script_ver="2.0.2"

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

    #mv Aligned.toTranscriptome.out.bam ${reads1_fn}-${reads2_fn}_star_anno.bam
    set -x
    mv Log.final.out ${reads1_fn}-${reads2_fn}_star_Log.final.out
    set +x

    echo "* Upload results..."
    star_genome_bam=$(dx upload ${reads1_fn}-${reads2_fn}_star_genome.bam --details="{ $meta }" --property SW="$versions" --brief)
    star_anno_bam=$(dx upload ${reads1_fn}-${reads2_fn}_star_anno.bam     --details="{ $meta }" --property SW="$versions" --brief)
    star_log=$(dx upload ${reads1_fn}-${reads2_fn}_star_Log.final.out --property SW="$versions" --brief)

    dx-jobutil-add-output star_genome_bam "$star_genome_bam" --class=file
    dx-jobutil-add-output star_anno_bam "$star_anno_bam" --class=file
    dx-jobutil-add-output star_log "$star_log" --class=file
    dx-jobutil-add-output metadata "$meta" --class=string

    echo "* Finished."
}
