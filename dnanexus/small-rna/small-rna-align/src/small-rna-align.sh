#!/bin/bash
# small-rna-align.sh

script_name="small-rna-align.sh"
script_ver="1.0.2"

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
    set -x
    libraryComment="@CO\tLIBID:${library_id}"
    echo -e ${libraryComment} > COfile.txt
    cat out/*_bamCommentLines.txt >> COfile.txt
    echo `cat COfile.txt`
    set +x

    echo "* Map reads..."
    set -x
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
    set +x

    echo "* Prepare metadata..."
    meta=''
    if [ -f /usr/bin/qc_metrics.py ]; then
        meta=`qc_metrics.py -n STAR_log_final -f Log.final.out`
    fi

    #mv Aligned.toTranscriptome.out.bam ${reads1_fn}-${reads2_fn}_star_anno.bam
    set -x
    mv Log.final.out ${reads_fn}_star_Log.final.out
    set +x

    echo "* Upload results..."
    genome_bam=$(dx upload ${reads_fn}_star_genome.bam --details "{ $meta }" --property SW="$versions" --brief)
    star_log=$(dx upload ${reads_fn}_star_Log.final.out --brief)

    dx-jobutil-add-output genome_bam "$genome_bam" --class=file
    dx-jobutil-add-output star_log "$star_log" --class=file
    dx-jobutil-add-output metadata "{ $meta }" --class=string

    echo "* Finished."
}
