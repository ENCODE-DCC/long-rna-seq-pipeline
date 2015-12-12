#!/bin/bash
# small-rna-align.sh

main() {
    # Now in resources/usr/bin
    #echo "* Download and install STAR..."
    #git clone https://github.com/alexdobin/STAR
    #(cd STAR; git checkout tags/STAR_2.4.2a)
    #(cd STAR; make)
    #wget https://github.com/samtools/samtools/archive/0.1.19.tar.gz

    # If available, will print tool versions to stderr and json string to stdout
    versions=''
    if [ -f /usr/bin/tool_versions.py ]; then 
        versions=`tool_versions.py --dxjson dnanexus-executable.json`
    fi

    echo "* Value of reads: '$reads'"
    echo "* Value of star_index: '$star_index'"
    echo "* Value of library_id: '$library_id'"
    echo "* Value of nthreads: '$nthreads'"

    #echo "* Download files..."
    outfile_name=""
    concat=""
    rm -f concat.fq
    for ix in ${!reads[@]}
    do
        file_root=`dx describe "${reads[$ix]}" --name`
        file_root=${file_root%.fastq.gz}
        file_root=${file_root%.fq.gz}
        if [ "${outfile_name}" == "" ]; then
            outfile_name="${file_root}"
        else
            outfile_name="${file_root}_${outfile_name}"
            if [ "${concat}" == "" ]; then
                outfile_name="${outfile_name}_concat" 
                concat="s concatenated as"
            fi
        fi
        echo "* Downloading and concatenating ${file_root}.fq.gz file..."
        dx download "${reads[$ix]}" -o - | gunzip >> concat.fq
    done
    mv concat.fq ${outfile_name}.fq
    echo "* Gzipping file..."
    gzip ${outfile_name}.fq
    echo "* Fastq${concat} file: '${outfile_name}.fq.gz'"
    reads_root=${outfile_name}
    ls -l ${reads_root}.fq.gz
    bam_root="${reads_root}_srna_star"
    if [ -f /usr/bin/parse_property.py ]; then
        new_root=`parse_property.py --job "${DX_JOB_ID}" --root_name --quiet`
        if [ "$new_root" != "" ]; then
            bam_root="${new_root}_srna_star"
        fi
    fi
    echo "* Alignments file will be: '${bam_root}.bam'"

    echo "* Downloading and extracting star index archive..."
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
    STAR --genomeDir out --readFilesIn ${reads_root}.fq.gz --readFilesCommand zcat              \
        --runThreadN ${nthreads} --outFilterMultimapNmax 20 --alignIntronMax 1                  \
        --clip3pAdapterSeq TGGAATTCTC --clip3pAdapterMMp 0.1 --outFilterMismatchNoverLmax 0.03  \
        --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 16  \
        --outSAMheaderCommentFile COfile.txt --outSAMheaderHD @HD VN:1.4 SO:coordinate          \
        --genomeLoad NoSharedMemory --outSAMunmapped Within --outSAMtype BAM SortedByCoordinate \
        --quantMode GeneCounts --alignSJDBoverhangMin 1000 --limitBAMsortRAM 60000000000
        
    mv Aligned.sortedByCoord.out.bam ${bam_root}.bam
    mv ReadsPerGene.out.tab ${bam_root}_quant.tsv
    mv Log.final.out ${bam_root}_Log.final.out
    set +x
    ls -l ${bam_root}.bam
    ls -l ${bam_root}_quant.tsv

    echo "* Prepare metadata..."
    meta=''
    if [ -f /usr/bin/qc_metrics.py ]; then
        meta=`qc_metrics.py -n STAR_log_final -f ${bam_root}_Log.final.out`
    fi

    echo "* Upload results..."
    srna_bam=$(dx upload ${bam_root}.bam --details "{ $meta }" --property SW="$versions" --brief)
    srna_quant=$(dx upload ${bam_root}_quant.tsv --details "{ $meta }" --property SW="$versions" --brief)
    star_log=$(dx upload ${bam_root}_Log.final.out --brief)

    dx-jobutil-add-output srna_bam "$srna_bam" --class=file
    dx-jobutil-add-output srna_quant "$srna_quant" --class=file
    dx-jobutil-add-output star_log "$star_log" --class=file
    dx-jobutil-add-output metadata "{ $meta }" --class=string

    echo "* Finished."
}
