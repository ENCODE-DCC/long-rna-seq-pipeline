#!/bin/bash
# prep-rsem 1.0.1

main() {
    # Now in resources/usr/bin
    #echo "* Dowload and install RSEM..."
    #git clone https://github.com/bli25wisc/RSEM.git
    # (cd RSEM; git checkout tags/v1.2.19)
    ### get correct commit bit from submodule
    #(cd RSEM; make)

    echo "*****"
    echo "* Running: prep-rsem.sh [v1.0.1]"
    echo "* RSEM version: "`rsem-calculate-expression --version | awk '{print $5}'`
    echo "*****"

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
    mkdir out
    rsem-prepare-reference --no-polyA --gtf ${annotation_fn}.gtf ${ref} out/rsem

    # Attempt to make bamCommentLines.txt, which should be reviewed. NOTE tabs handled by assignment.
    echo "* Create bam header..."
    refComment="@CO\tREFID:$(basename ${genome_fn})"
    annotationComment="@CO\tANNID:$(basename ${annotation_fn})"
    spikeInComment="@CO\tSPIKEID:${spike_in_fn}"
    echo -e ${refComment} > out/rsem_bamCommentLines.txt
    echo -e ${annotationComment} >> out/rsem_bamCommentLines.txt
    echo -e ${spikeInComment} >> out/rsem_bamCommentLines.txt

    echo `cat "out/rsem_bamCommentLines.txt"`

    echo "* Tar and upload..."
    echo `ls out/*`
    tar -czf ${genome_fn}_${annotation_fn}_rsemIndex.tgz out/*

    rsem_index=$(dx upload ${genome_fn}_${annotation_fn}_rsemIndex.tgz --brief)
    dx-jobutil-add-output rsem_index $rsem_index --class=file
    echo "* Finished."
}
