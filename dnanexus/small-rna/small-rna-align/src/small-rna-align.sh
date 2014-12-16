#!/bin/bash
# small-rna-align 0.0.1

main() {
    # Now in resources/usr/bin
    #echo "* Download and install STAR..."
    #git clone https://github.com/alexdobin/STAR
    #(cd STAR; git checkout tags/STAR_2.4.0g1)
    #(cd STAR; make)
    #wget https://github.com/samtools/samtools/archive/0.1.19.tar.gz

    echo "*****"
    echo "* Running: small-rna-align.sh [v0.0.1]"
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
    STAR --genomeDir out --readFilesIn ${reads_fn}.fastq.gz --readFilesCommand zcat \
        --runThreadN ${nthreads} --outFilterMultimapNmax 20 --alignIntronMax 1 \
        --clip3pAdapterSeq TGGAATTCTC --clip3pAdapterMMp 0.1 --outFilterMismatchNoverLmax 0.05 \
        --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 16 \
        --outSAMheaderCommentFile COfile.txt --outSAMheaderHD @HD VN:1.4 SO:coordinate   \
        --genomeLoad NoSharedMemory --outSAMunmapped Within --outSAMtype BAM SortedByCoordinate
        
    #echo "* Index genome bam..."
    # Note: no longer making unused index
    mv Aligned.sortedByCoord.out.bam ${reads_fn}_star_genome.bam
    #samtools index ${reads_fn}_star_genome.bam

    mv Log.final.out ${reads_fn}_star_Log.final.out

    echo "* Upload results..."
    star_log=$(dx upload ${reads_fn}_star_Log.final.out --brief)
    genome_bam=$(dx upload ${reads_fn}_star_genome.bam --brief)
    #genome_bai=$(dx upload ${reads_fn}_star_genome.bam.bai --brief)

    dx-jobutil-add-output star_log "$star_log" --class=file
    dx-jobutil-add-output genome_bam "$genome_bam" --class=file
    #dx-jobutil-add-output genome_bai "$genome_bai" --class=file
    echo "* Finished."
}
