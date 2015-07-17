#!/bin/bash
# prep-star.sh

main() {
    # Now in resources/usr/bin
    #echo "* Download and install STAR..."
    #git clone https://github.com/alexdobin/STAR
    #(cd STAR; git checkout tags/STAR_2.4.0h)
    #(cd STAR; make)

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

    echo "* Build index..."
    set -x
    mkdir out
    STAR --runMode genomeGenerate --genomeFastaFiles ${genome_fn}.fa ${spike_in_fn}.fa \
         --sjdbOverhang 100 --sjdbGTFfile ${annotation_fn}.gtf --runThreadN 6 --genomeDir out/ \
         --outFileNamePrefix out
    set +x

    # Attempt to make bamCommentLines.txt, which should be reviewed. NOTE tabs handled by assignment.
    echo "* Create bam header..."
    set -x
    refComment="@CO\tREFID:$(basename ${genome_fn})"
    annotationComment="@CO\tANNID:$(basename ${annotation_fn})"
    spikeInComment="@CO\tSPIKEID:${spike_in_fn}"
    echo -e ${refComment} > out/star_bamCommentLines.txt
    echo -e ${annotationComment} >> out/star_bamCommentLines.txt
    echo -e ${spikeInComment} >> out/star_bamCommentLines.txt

    echo `cat "out/star_bamCommentLines.txt"`
    set +x

    echo "* Tar and upload results..."
    echo "/"
    echo `ls`
    echo "out"
    echo `ls out/`
    set -x
    tar -czf ${genome_fn}_${annotation_fn}_starIndex.tgz out/
    set +x

    star_index=$(dx upload ${genome_fn}_${annotation_fn}_starIndex.tgz --property SW="$versions" --brief)

    # The following line(s) use the utility dx-jobutil-add-output to format and
    # add output variables to your job's output as appropriate for the output
    # class.  Run "dx-jobutil-add-output -h" for more information on what it
    # does.

    dx-jobutil-add-output star_index $star_index --class=file
    echo "* Finished."
}
