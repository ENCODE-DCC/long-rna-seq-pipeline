#!/bin/bash
# concat-fastqs.sh

script_name="concat-fastqs.sh"
script_ver="1.0.4"

main() {
    # If available, will print tool versions to stderr and json string to stdout
    versions=''
    if [ -f /usr/bin/tool_versions.py ]; then 
        versions=`tool_versions.py --applet $script_name --appver $script_ver`
    fi

    echo "* Concat id: '${concat_id}'"
    #echo "* Value of trna_annoation: '$trna_annotation'"
    #echo "* Value of spike_in: '$spike_in'"

    #echo "* Download and concatentating files..."
    outfile_root="${concat_id}_concat"
    for ix in ${!reads_set[@]}
    do
        file_root=`dx describe "${reads_set[$ix]}" --name`
        file_root=${file_root%.fastq.gz}
        file_root=${file_root%.fq.gz}
        echo "* Download and concatentating '${file_root}.fq.gz'..."
        set -x
        dx download "${reads_set[$ix]}" -o - | gunzip >> concat.fq
        set +x
        outfile_root="${file_root}_${outfile_root}"
    done
    echo "* Expected name of output file: '${outfile_root}.fq.gz'"

    echo "* Gzipping file..."
    set -x
    mv >> concat.fq ${outfile_root}.fq
    gzip ${outfile_root}.fq
    set +x

    echo "* Upload results..."
    reads=$(dx upload ${outfile_root}.fq.gz --property SW="$versions" --brief)
    dx-jobutil-add-output reads "$reads" --class=file
    echo "* Finished."
}
