#!/bin/bash
# concat-fastqs 0.0.2

main() {
    echo "*****"
    echo "* Running: concat-fastqs"
    echo "*****"

    echo "* Expected name of ouutput file: '${outfile_root}.fastq.gz'"
    #echo "* Value of trna_annoation: '$trna_annotation'"
    #echo "* Value of spike_in: '$spike_in'"

    echo "* Download files..."
    for ix in ${!fastq_files[@]}
    do
        filename=`dx describe "${fastq_files[$ix]}" --name | cut -d'.' -f1`
        dx download "${fastq_files[$ix]}" -o - | gunzip > "$filename".fq
    done

    echo "* Concatenating files..."
    cat *.fq > ${outfile_root}.fq

    echo "* Gzipping file..."
    gzip ${outfile_root}.fq

    echo "* Upload results..."
    combined_fastq=$(dx upload ${outfile_root}.fq.gz --brief)
    dx-jobutil-add-output combined_fastq "$combined_fastq" --class=file
    echo "* Finished."
}
