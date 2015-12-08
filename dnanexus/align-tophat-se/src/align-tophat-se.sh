#!/bin/bash
# align-tophat-se.sh

main() {
    # Now in resources/usr/bin
    ## install tophat 2.0.8
    #wget http://ccb.jhu.edu/software/tophat/downloads/tophat-2.0.8.Linux_x86_64.tar.gz
    ## install bowtie2_2.1.0
    #wget http://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.1.0/bowtie2-2.1.0-linux-x86_64.zip/download
    #wget https://github.com/samtools/samtools/archive/0.1.19.tar.gz
    #wget wget https://github.com/xweigit/xweiEncodeScripts/archive/v1.0.tar.gz

    # If available, will print tool versions to stderr and json string to stdout
    versions=''
    if [ -f /usr/bin/tool_versions.py ]; then 
        versions=`tool_versions.py --dxjson dnanexus-executable.json`
    fi

    echo "* Value of reads: '$reads'"
    echo "* Value of tophat_index: '$tophat_index'"
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
    bam_root="${reads_root}"
    if [ -f /usr/bin/parse_property.py ]; then
        new_root=`parse_property.py -f "'${reads[0]}'" --project "${DX_PROJECT_CONTEXT_ID}" --root_name --quiet`
        if [ "$new_root" != "" ]; then
            bam_root="${new_root}"
        fi
    fi

    echo "* Downloading and extracting TopHat index archive..."
    dx download "$tophat_index" -o tophat_index.tgz
  
    echo "* ===== Calling DNAnexus and ENCODE independent script... ====="
    set -x
    lrna-align-tophat-se.sh tophat_index.tgz ${reads_root}.fq.gz "$library_id" $nthreads $bam_root
    set +x
    echo "* ===== Returned from dnanexus and encodeD independent script ====="
    bam_root="${bam_root}_tophat"

    echo "* Prepare metadata..."
    qc_stats=''
    reads=0
    if [ -f /usr/bin/qc_metrics.py ]; then
        qc_stats=`qc_metrics.py -n samtools_flagstats -f ${bam_root}_flagstat.txt`
        reads=`qc_metrics.py -n samtools_flagstats -f ${bam_root}_flagstat.txt -k total`
    fi

    echo "* Upload results..."
    tophat_bam=$(dx upload ${bam_root}.bam --details="{ $qc_stats }" --property reads="$reads" \
                                                                     --property SW="$versions" --brief)
    tophat_flagstat=$(dx upload ${bam_root}_flagstat.txt --details="{ $qc_stats }" --property reads="$reads" \
                                                                                   --property SW="$versions" --brief)

    dx-jobutil-add-output tophat_bam "$tophat_bam" --class=file
    dx-jobutil-add-output tophat_flagstat "$tophat_flagstat" --class=file

    dx-jobutil-add-output reads "$reads" --class=string
    dx-jobutil-add-output metadata "{ $qc_stats }" --class=string

    echo "* Finished."
}
