#!/bin/bash
# small-rna-prep-star 0.0.1

main() {
    # Now in resources/usr/bin
    #echo "* Download and install STAR..."
    #git clone https://github.com/alexdobin/STAR
    #(cd STAR; git checkout tags/STAR_2.4.0g1)
    #(cd STAR; make)

    echo "*****"
    echo "* Running: small-rna-prep-star.sh [v0.0.1]"
    echo "* STAR version: "`STAR --version | awk '{print $1}' | cut -d _ -f 2-`
    echo "*****"

    echo "* Value of genome: '$genome'"

    # The following line(s) use the dx command-line tool to download your file
    # inputs to the local file system using variable names for the filenames. To
    # recover the original filenames, you can use the output of "dx describe
    # "$variable" --name".

    echo "* Download files..."
    genome_fn=`dx describe "$genome" --name`
    genome_fn=${genome_fn%.fasta.gz}
    genome_fn=${genome_fn%.fa.gz}
    dx download "$genome" -o "$genome_fn".fa.gz
    gunzip "$genome_fn".fa.gz
    ref="$genome_fn".fa
    echo "* Reference file(s): '$ref'"

    # Fill in your application code here.

    echo "* Build index..."
    mkdir out
    STAR --runMode genomeGenerate --genomeFastaFiles ${genome_fn}.fa --runThreadN 8 \
         --genomeDir out/ --outFileNamePrefix out

    # Attempt to make bamCommentLines.txt, which should be reviewed. NOTE tabs handled by assignment.
    echo "* Create bam header..."
    refComment="@CO\tREFID:$(basename ${genome_fn})"
    echo -e ${refComment} > out/star_bamCommentLines.txt
    echo `cat "out/star_bamCommentLines.txt"`

    echo "* Tar and upload results..."
    echo "out"
    echo `ls out/`
    tar -czf ${genome_fn}_starIndex.tgz out/

    star_index=$(dx upload ${genome_fn}_starIndex.tgz --brief)

    # The following line(s) use the utility dx-jobutil-add-output to format and
    # add output variables to your job's output as appropriate for the output
    # class.  Run "dx-jobutil-add-output -h" for more information on what it
    # does.

    dx-jobutil-add-output star_index $star_index --class=file
    echo "* Finished."
}
