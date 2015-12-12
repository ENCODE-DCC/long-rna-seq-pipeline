#!/bin/bash
# rampage-align-pe.sh

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
        versions=`tool_versions.py --dxjson dnanexus-executable.json`
    fi

    echo "* Value of reads1: '$reads1'"
    echo "* Value of reads2: '$reads2'"
    echo "* Value of star_index: '$star_index'"
    echo "* Value of library_id: '$library_id'"
    echo "* Number of threads (default 8): '$nthreads'"

    #echo "* Download files..."
    outfile_name=""
    concat=""
    rm -f concat.fq
    for ix in ${!reads1[@]}
    do
        file_root=`dx describe "${reads1[$ix]}" --name`
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
        echo "* Downloading concatenating ${file_root}.fq.gz file..."
        dx download "${reads1[$ix]}" -o - | gunzip >> concat.fq
    done
    mv concat.fq ${outfile_name}.fq
    echo "* Gzipping file..."
    gzip ${outfile_name}.fq
    echo "* Reads1 fastq${concat} file: '${outfile_name}.fq.gz'"
    reads1_root=${outfile_name}
    ls -l ${reads1_root}.fq.gz

    outfile_name=""
    concat=""
    rm -f concat.fq
    for ix in ${!reads2[@]}
    do
        file_root=`dx describe "${reads2[$ix]}" --name`
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
        dx download "${reads2[$ix]}" -o - | gunzip >> concat.fq
    done
    mv concat.fq ${outfile_name}.fq
    echo "* Gzipping file..."
    gzip ${outfile_name}.fq
    echo "* Reads2 fastq${concat} file: '${outfile_name}.fq.gz'"
    ls -l ${outfile_name}.fq.gz
    reads2_root=${outfile_name}
    ls -l ${reads2_root}.fq.gz
    bam_root="${reads1_root}_${reads2_root}"
    if [ -f /usr/bin/parse_property.py ]; then
        new_root=`parse_property.py --job "${DX_JOB_ID}" --root_name --quiet`
        if [ "$new_root" != "" ]; then
            bam_root="${new_root}"
        fi
    fi

    echo "* Downloading star index archive..."
    dx download "$star_index" -o star_index.tgz
    
    # DX/ENCODE independent script is found in resources/usr/bin
    echo "* ===== Calling DNAnexus and ENCODE independent script... ====="
    set -x
    ram-align-star-pe.sh star_index.tgz ${reads1_root}.fq.gz ${reads2_root}.fq.gz "$library_id" $nthreads $bam_root
    set +x
    echo "* ===== Returned from dnanexus and encodeD independent script ====="
    bam_root="${bam_root}_rampage_star"

    echo "* Prepare metadata..."
    qc_stats=''
    reads=0
    if [ -f /usr/bin/qc_metrics.py ]; then
        qc_stats=`qc_metrics.py -n STAR_log_final -f ${bam_root}_Log.final.out`
        meta=`qc_metrics.py -n samtools_flagstats -f ${bam_root}_marked_flagstat.txt`
        reads=`qc_metrics.py -n samtools_flagstats -f ${bam_root}_marked_flagstat.txt -k total`
        qc_stats=`echo $qc_stats, $meta`
    fi

    echo "* Upload results..."
    rampage_marked_bam=$(dx upload ${bam_root}_marked.bam --details "{ $qc_stats }" --property reads="$reads" \
                                                                                    --property SW="$versions" --brief)
    rampage_star_log=$(dx upload ${bam_root}_Log.final.out --details "{ $qc_stats }" --property SW="$versions" --brief)
    rampage_flagstat=$(dx upload ${bam_root}_marked_flagstat.txt --details="{ $qc_stats }" --property reads="$reads" \
                                                                                           --property SW="$versions" --brief)

    dx-jobutil-add-output rampage_marked_bam "$rampage_marked_bam" --class=file
    dx-jobutil-add-output rampage_star_log "$rampage_star_log" --class=file
    dx-jobutil-add-output rampage_flagstat "$rampage_flagstat" --class=file

    dx-jobutil-add-output reads "$reads" --class=string
    dx-jobutil-add-output metadata "{ $qc_stats }" --class=string
    echo "* Finished."
}
