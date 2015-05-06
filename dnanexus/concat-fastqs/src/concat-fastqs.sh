#!/bin/bash
# concat-fastqs.sh

script_name="concat-fastqs.sh"
script_ver="1.0.3"

main() {
    # If available, will print tool versions to stderr and json string to stdout
    versions=''
    if [ -f /usr/bin/tool_versions.py ]; then 
        versions=`tool_versions.py --applet $script_name --appver $script_ver`
    fi
    #echo "*****"
    #echo "* Running: concat-fastqs [v1.0.2]"
    #echo "*****"

    echo "* Concat id: '${concat_id}'"
    #echo "* Value of trna_annoation: '$trna_annotation'"
    #echo "* Value of spike_in: '$spike_in'"

    echo "* Download files..."
    outfile_root="${concat_id}_concat"
    for ix in ${!reads_set[@]}
    do
        filename=`dx describe "${reads_set[$ix]}" --name | cut -d'.' -f1`
        dx download "${reads_set[$ix]}" -o - | gunzip > ${filename}.fastq
        file_root=${filename%.fastq.gz}
        file_root=${filename%.fq.gz}
        outfile_root="${file_root}_${outfile_root}"
    done
    echo "* Expected name of output file: '${outfile_root}.fq.gz'"

    echo "* Concatenating files..."
    # Note: must maintain order so that paired reads maintain order.
    set -x
    rm -f ${outfile_root}.fq
    for ix in ${!reads_set[@]}
    do
        filename=`dx describe "${reads_set[$ix]}" --name | cut -d'.' -f1`
        cat ${filename}.fastq >> ${outfile_root}.fq
    done
    set +x
    
    echo "* Gzipping file..."
    set -x
    gzip ${outfile_root}.fq
    set +x

    echo "* Upload results..."
    reads=$(dx upload ${outfile_root}.fq.gz --property SW="$versions" --brief)
    dx-jobutil-add-output reads "$reads" --class=file
    echo "* Finished."
}
