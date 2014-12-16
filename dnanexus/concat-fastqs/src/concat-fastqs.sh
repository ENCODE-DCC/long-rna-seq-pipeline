#!/bin/bash
# concat-fastqs 1.0.1

main() {
    echo "*****"
    echo "* Running: concat-fastqs [v1.0.1]"
    echo "*****"

    echo "* Expected name of output file: '${outfile_root}.fastq.gz'"
    #echo "* Value of trna_annoation: '$trna_annotation'"
    #echo "* Value of spike_in: '$spike_in'"

    echo "* Download files..."
    for ix in ${!reads_set[@]}
    do
        filename=`dx describe "${reads_set[$ix]}" --name | cut -d'.' -f1`
        dx download "${reads_set[$ix]}" -o - | gunzip > ${filename}.fastq
    done

    echo "* Concatenating files..."
    # Note: must maintain order so that paired reads maintain order.
    rm -f ${outfile_root}.fq
    for ix in ${!reads_set[@]}
    do
        filename=`dx describe "${reads_set[$ix]}" --name | cut -d'.' -f1`
        cat ${filename}.fastq >> ${outfile_root}.fq
    done
    
    echo "* Gzipping file..."
    gzip ${outfile_root}.fq

    echo "* Upload results..."
    reads=$(dx upload ${outfile_root}.fq.gz --brief)
    dx-jobutil-add-output reads "$reads" --class=file
    echo "* Finished."
}
