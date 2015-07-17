#!/bin/bash
# prep-rsem.sh

main() {
    # Now in resources/usr/bin
    #echo "* Dowload and install RSEM..."
    #git clone https://github.com/bli25wisc/RSEM.git
    # (cd RSEM; git checkout tags/v1.2.19)
    ### get correct commit bit from submodule
    #(cd RSEM; make)

    # If available, will print tool versions to stderr and json string to stdout
    versions=''
    if [ -f /usr/bin/tool_versions.py ]; then 
        versions=`tool_versions.py --dxjson dnanexus-executable.json`
    fi

    echo "* Value of annotations: '$annotations'"
    echo "* Value of genome: '$genome'"
    echo "* Value of spike_in: '$spike_in'"

    # The following line(s) use the dx command-line tool to download your file
    # inputs to the local file system using variable names for the filenames. To
    # recover the original filenames, you can use the output of "dx describe
    # "$variable" --name".

    echo "* Download files..."
    annotation_fn=`dx describe "$annotations" --name`
    annotation_fn=${annotation_fn%.gtf.gz}
    dx download "$annotations" -o "$annotation_fn".gtf.gz
    gunzip "$annotation_fn".gtf.gz

    genome_fn=`dx describe "$genome" --name`
    genome_fn=${genome_fn%.fasta.gz}
    genome_fn=${genome_fn%.fa.gz}
    dx download "$genome" -o "$genome_fn".fa.gz
    gunzip "$genome_fn".fa.gz
    ref="$genome_fn".fa

    if [ -n "$spike_in" ]
    then
        spike_in_fn=`dx describe "$spike_in" --name`
        spike_in_fn=${spike_in_fn%.fasta.gz}
        spike_in_fn=${spike_in_fn%.fa.gz}
        dx download "$spike_in" -o "$spike_in_fn".fa.gz
        gunzip "$spike_in_fn".fa.gz
        ref="${ref},${spike_in_fn}.fa"
    fi
    echo "* Reference file(s): '$ref'"

    # Fill in your application code here.

    echo "* Prepare reference..."
    set -x
    mkdir out
    rsem-prepare-reference --no-polyA --gtf ${annotation_fn}.gtf ${ref} out/rsem
    set +x

    # Attempt to make bamCommentLines.txt, which should be reviewed. NOTE tabs handled by assignment.
    echo "* Create bam header..."
    set -x
    refComment="@CO\tREFID:$(basename ${genome_fn})"
    annotationComment="@CO\tANNID:$(basename ${annotation_fn})"
    spikeInComment="@CO\tSPIKEID:${spike_in_fn}"
    echo -e ${refComment} > out/rsem_bamCommentLines.txt
    echo -e ${annotationComment} >> out/rsem_bamCommentLines.txt
    echo -e ${spikeInComment} >> out/rsem_bamCommentLines.txt

    echo `cat "out/rsem_bamCommentLines.txt"`
    set +x

    echo "* Tar and upload..."
    echo `ls out/*`
    set -x
    tar -czf ${genome_fn}_${annotation_fn}_rsemIndex.tgz out/*
    set +x

    rsem_index=$(dx upload ${genome_fn}_${annotation_fn}_rsemIndex.tgz --property SW="$versions" --brief)
    dx-jobutil-add-output rsem_index $rsem_index --class=file
    echo "* Finished."
}
