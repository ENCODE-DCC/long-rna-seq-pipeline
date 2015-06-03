#!/bin/bash
# rampage-align-pe.sh

script_name="rampage-align-pe.sh"
script_ver="1.0.1"

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
    
    aligned_root=${reads1_fn}-${reads2_fn}_rampage_star
    echo "* Rampage Aligned file root:'${aligned_root}'"

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
    STAR --genomeDir out --readFilesIn ${reads1_fn}.fastq.gz ${reads2_fn}.fastq.gz       \
        --readFilesCommand zcat --runThreadN ${nthreads} --genomeLoad NoSharedMemory      \
        --outFilterMultimapNmax 500 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1        \
        --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04                   \
        --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000              \
        --outSAMheaderCommentFile COfile.txt --outSAMheaderHD @HD VN:1.4 SO:coordinate        \
        --outSAMunmapped Within --outFilterType BySJout --outSAMattributes NH HI AS NM MD      \
        --outFilterScoreMinOverLread 0.85 --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
        --clip5pNbases 6 15 --seedSearchStartLmax 30 --outSAMtype BAM SortedByCoordinate         \
        --limitBAMsortRAM 60000000000
    set +x

    echo "* Marking PCR duplicates..."
    set -x
    STAR --inputBAMfile Aligned.sortedByCoord.out.bam --bamRemoveDuplicatesType UniqueIdentical \
        --runMode inputAlignmentsFromBAM --bamRemoveDuplicatesMate2basesN 15 \
        --outFileNamePrefix markdup. --limitBAMsortRAM 60000000000
    set +x

    echo "* Prepare metadata..."
    meta=''
    if [ -f /usr/bin/qc_metrics.py ]; then
        meta=`qc_metrics.py -n STAR_log_final -f Log.final.out`
    fi

    #mv Aligned.toTranscriptome.out.bam ${reads_fn}_star_anno.bam
    set -x
    mv Log.final.out ${aligned_root}_Log.final.out
    set +x

    echo "* Upload results..."
    set -x
    mv markdup.Processed.out.bam ${aligned_root}_marked.bam
    set +x

    # NOTE: adding meta 'details' ensures json is valid.  But details are not updatable so rely on QC property
    rampage_marked_bam=$(dx upload ${aligned_root}_marked.bam --details "{ $meta }" --property QC="{ $meta }" --property SW="$versions" --brief)
    rampage_star_log=$(dx upload ${aligned_root}_Log.final.out --property SW="$versions" --brief)

    dx-jobutil-add-output rampage_marked_bam "$rampage_marked_bam" --class=file
    dx-jobutil-add-output rampage_star_log "$rampage_star_log" --class=file
    dx-jobutil-add-output metadata "{ $meta }" --class=string
    echo "* Finished."
}
