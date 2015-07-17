#!/bin/bash
# merge-annotation.sh

main() {
    # Now in resources/usr/bin
    #wget https://github.com/ENCODE-DCC/long-rna-seq-pipeline/tree/master/DAC/GTF.awk

    # If available, will print tool versions to stderr and json string to stdout
    versions=''
    if [ -f /usr/bin/tool_versions.py ]; then 
        versions=`tool_versions.py --dxjson dnanexus-executable.json`
    fi

    echo "* Value of gene_annotation: '$gene_annotation'"
    echo "* Value of trna_annoation: '$trna_annotation'"
    echo "* Value of spike_in: '$spike_in'"

    echo "* Download files..."
    gene_fn=`dx describe "$gene_annotation" --name`
    gene_fn=${gene_fn%.gtf.gz}
    dx download "$gene_annotation" -o "$gene_fn".gtf.gz
    gunzip "$gene_fn".gtf.gz

    trna_fn=`dx describe "$trna_annotation" --name`
    trna_fn=${trna_fn%.gtf.gz}
    dx download "$trna_annotation" -o "$trna_fn".gtf.gz
    gunzip "$trna_fn".gtf.gz

    spike_in_fn=`dx describe "$spike_in" --name`
    spike_in_fn=${spike_in_fn%.fasta.gz}
    spike_in_fn=${spike_in_fn%.fa.gz}
    dx download "$spike_in" -o "$spike_in_fn".fa.gz
    gunzip "$spike_in_fn".fa.gz

    # Fill in your application code here.

    echo "* Merge GTF files..."
    set -x
    out_fn="$gene_fn"-"tRNAs"-"$spike_in_fn".gtf
    awk -f /usr/bin/GTF.awk ${gene_fn}.gtf ${trna_fn}.gtf ${spike_in_fn}.fa > ${out_fn}
    gzip ${out_fn}
    set +x

    echo "* Upload results..."
    combined_gtf=$(dx upload $out_fn.gz --property SW="$versions" --brief)
    dx-jobutil-add-output combined_gtf "$combined_gtf" --class=file
    echo "* Finished."
}
