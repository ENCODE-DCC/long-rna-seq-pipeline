#!/bin/bash
# prep-tophat.sh

script_name="prep-tophat.sh"
script_ver="1.0.1"

main() {
    # Now in resources/usr/bin
    ## install tophat 2.0.8
    #wget http://ccb.jhu.edu/software/tophat/downloads/tophat-2.0.8.Linux_x86_64.tar.gz
    ## install bowtie2_2.1.0
    #wget http://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.1.0/bowtie2-2.1.0-linux-x86_64.zip/download

    # If available, will print tool versions to stderr and json string to stdout
    versions=''
    if [ -f /usr/bin/tool_versions.py ]; then 
        versions=`tool_versions.py --applet $script_name --appver $script_ver`
    fi
    #echo "*****"
    #echo "* Running: prep_tophat.sh [v1.0.0]"
    #echo "* TopHat version: "`tophat -v | awk '{print $2}'`
    #echo "* bowtie2 version: "`bowtie2 --version 2>&1 | grep bowtie | awk '{print $3}'`
    #echo "*****"

    echo "* Value of annotations: '$annotations'"
    echo "* Value of genome: '$genome'"
    echo "* Value of spike_in: '$spike_in'"

    echo "* Download files..."
    annotation_fn=`dx describe "$annotations" --name`
    annotation_fn=${annotation_fn%.gtf.gz}
    dx download "$annotations" -o "$annotation_fn".gtf.gz
    gunzip "$annotation_fn".gtf.gz
    anno_prefix=${annotation_fn}

    genome_fn=`dx describe "$genome" --name`
    genome_fn=${genome_fn%.fasta.gz}
    genome_fn=${genome_fn%.fa.gz}
    dx download "$genome" -o "$genome_fn".fa.gz
    gunzip "$genome_fn".fa.gz
    ref="$genome_fn".fa
    geno_prefix=${genome_fn}

    dx download file-BKyfb080ZZ0P4jQFVGB01966 -o tiny.fq.gz
    gunzip tiny.fq.gz

    if [ -n "$spike_in" ]; then
        spike_in_fn=`dx describe "$spike_in" --name`
        spike_in_fn=${spike_in_fn%.fasta.gz}
        spike_in_fn=${spike_in_fn%.fa.gz}
        dx download "$spike_in" -o "$spike_in_fn".fa.gz
        gunzip "$spike_in_fn".fa.gz
        ref="${ref},${spike_in_fn}.fa"
        geno_prefix="${genome_fn}-${spike_in_fn}"
    fi
    echo "* Reference file(s): '$ref'"
    echo "* Value of geno_prefix: '$geno_prefix'"
    echo "* Value of anno_prefix: '$anno_prefix'"

    # Fill in your application code here.

    set -x
    mkdir out
    #bowtie2-build --offrate 3 -f ${ref} out/$geno_prefix  ### Definitely has an effect on results
    bowtie2-build -f ${ref} out/$geno_prefix
    # make sure the combined fa file is preserved in the archive, so that it isn't rebuilt each time
    bowtie2-inspect out/$geno_prefix > out/$geno_prefix.fa
    set +x

    # Attempt to make bamCommentLines.txt. NOTE tabs handled by assignment.
    echo "* Create bam header..."
    set -x
    refComment="@CO\tREFID:${genome_fn}"
    annotationComment="@CO\tANNID:${annotation_fn}"
    echo -e ${refComment} > out/${geno_prefix}_bamCommentLines.txt
    echo -e ${annotationComment} >> out/${geno_prefix}_bamCommentLines.txt
    if [ -n "$spike_in" ]; then
        spikeInComment="@CO\tSPIKEID:${spike_in_fn}"
        echo -e ${spikeInComment} >> out/${geno_prefix}_bamCommentLines.txt
    fi

    echo `cat "out/${geno_prefix}_bamCommentLines.txt"`
    set +x

    echo "* Run a 'quicky' tophat to generate index..."
    set -x
    tophat --no-discordant --no-mixed -p 8 -z0 --min-intron-length 20 --max-intron-length 1000000 \
           --read-mismatches 4 --read-edit-dist 4 --max-multihits 20 --library-type fr-firststrand \
           --GTF "$annotation_fn".gtf --no-coverage-search \
           --transcriptome-index=out/${anno_prefix} out/${geno_prefix} tiny.fq
    set +x

    echo "* Tar and upload results..."
    echo `ls out/`
    set -x
    tar -czf ${genome_fn}_${annotation_fn}_tophatIndex.tgz out/${geno_prefix}* out/${anno_prefix}*
    set +x

    tophat_index=$(dx upload ${genome_fn}_${annotation_fn}_tophatIndex.tgz --property SW="$versions" --brief)

    dx-jobutil-add-output tophat_index $tophat_index --class=file
    echo "* Finished."
}

