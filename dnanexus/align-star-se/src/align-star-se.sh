#!/bin/bash
# align-star-se 1.0.1

main() {
    # Now in resources/usr/bin
    #echo "* Download and install STAR..."
    #git clone https://github.com/alexdobin/STAR
    #(cd STAR; git checkout tags/STAR_2.4.0h)
    #(cd STAR; make)
    #wget https://github.com/samtools/samtools/archive/0.1.19.tar.gz

    echo "*****"
    echo "* Running: align-star-se.sh [v1.0.1]"
    echo "* STAR version: "`STAR --version | awk '{print $1}' | cut -d _ -f 2-`
    echo "* samtools version: "`samtools 2>&1 | grep Version | awk '{print $2}'`
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
    STAR --genomeDir out --readFilesIn ${reads_fn}.fastq.gz                         \
        --readFilesCommand zcat --runThreadN ${nthreads} --genomeLoad NoSharedMemory \
        --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1    \
        --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04              \
        --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000         \
        --outSAMheaderCommentFile COfile.txt --outSAMheaderHD @HD VN:1.4 SO:coordinate   \
        --outSAMunmapped Within --outFilterType BySJout --outSAMattributes NH HI AS NM MD \
        --outSAMstrandField intronMotif --outSAMtype BAM SortedByCoordinate                \
        --quantMode TranscriptomeSAM --sjdbScore 1

    #echo "* Index genome bam..."
    # Note: no longer making unused index
    mv Aligned.sortedByCoord.out.bam ${reads_fn}_star_genome.bam
    #samtools index ${reads_fn}_star_genome.bam

    echo "* Sorting annotation bam..."
    cat <( samtools view -H Aligned.toTranscriptome.out.bam ) \
        <( samtools view -@ ${nthreads} Aligned.toTranscriptome.out.bam | sort -S 60G -T ./ ) | \
        samtools view -@ ${nthreads} -bS - > ${reads_fn}_star_anno.bam

    #mv Aligned.toTranscriptome.out.bam ${reads_fn}_star_anno.bam
    mv Log.final.out ${reads_fn}_star_Log.final.out
    mv Log.out ${reads_fn}_star_Log.out

    echo "* Upload results..."
    #detail_log=$(dx upload ${reads_fn}_star_Log.out --brief)
    star_log=$(dx upload ${reads_fn}_star_Log.final.out --brief)
    star_genome_bam=$(dx upload ${reads_fn}_star_genome.bam --brief)
    #star_genome_bai=$(dx upload ${reads_fn}_star_genome.bam.bai --brief)
    star_anno_bam=$(dx upload ${reads_fn}_star_anno.bam --brief)

    #dx-jobutil-add-output detail_log "$detail_log" --class=file
    dx-jobutil-add-output star_log "$star_log" --class=file
    dx-jobutil-add-output star_genome_bam "$star_genome_bam" --class=file
    #dx-jobutil-add-output star_genome_bai "$star_genome_bai" --class=file
    dx-jobutil-add-output star_anno_bam "$star_anno_bam" --class=file
    echo "* Finished."
}
