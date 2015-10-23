#!/bin/bash
# quant-rsem.sh 

main() {
    # Now in resources/usr/bin
    # echo "* Dowload and install RSEM..."
    # git clone https://github.com/bli25wisc/RSEM.git    
    #wget https://github.com/deweylab/RSEM/archive/v1.2.23.tar.gz -O rsem.tgz >> install.log 2>&1
    #mkdir rsem
    #tar -xzf rsem.tgz -C rsem --strip-components=1
    #cd rsem
    #make >> install.log 2>&1
    #mv rsem-build-read-index /usr/bin/ >> install.log 2>&1
    #mv rsem-calculate-credibility-intervals /usr/bin/ >> install.log 2>&1
    #mv rsem-calculate-expression /usr/bin/ >> install.log 2>&1
    #mv rsem-parse-alignments /usr/bin/ >> install.log 2>&1
    #mv rsem-prepare-reference /usr/bin/ >> install.log 2>&1
    #mv rsem-run-em /usr/bin/ >> install.log 2>&1
    #mv rsem-run-gibbs /usr/bin/ >> install.log 2>&1
    #mv rsem_perl_utils.pm /usr/bin/ >> install.log 2>&1
    #cd .. 

    # If available, will print tool versions to stderr and json string to stdout
    versions=''
    if [ -f /usr/bin/tool_versions.py ]; then 
        versions=`tool_versions.py --dxjson dnanexus-executable.json`
    fi

    echo "* Value of annotation_bam: '$star_anno_bam'"
    echo "* Value of rsem_index: '$rsem_index'"
    echo "* Value of paired: '$paired_end'"
    #echo "* Value of stranded: '$stranded'"
    echo "* Random number seed: '$rnd_seed'"
    echo "* Value of nthreads: '$nthreads'"

    echo "* Download files..."
    bam_fn=`dx describe "$star_anno_bam" --name`
    bam_fn=${bam_fn%_star_anno.bam}
    bam_fn=${bam_fn%.bam}
    dx download "$star_anno_bam" -o ${bam_fn}.bam
    dx download "$rsem_index" -o rsem_index.tgz
    tar zxvf rsem_index.tgz

    # should be 'out/rsem'
    grp=`ls out/*.grp`
    index_prefix=${grp%.grp}
    echo "* Found index_prefix: '$index_prefix'"

    extraFlags=""
    if [ "$paired_end" == "true" ]; then
        echo '* Running for paired-end, stranded'
        extraFlags="--paired-end --forward-prob 0"
    else
        echo '* Running for unpaired, unstranded'
    fi

    #if [ "$stranded" == "true" ]
    #then
    #    echo '* Using stranded flag'
    #    extraFlags=${extraFlags}"--forward-prob 0"
    #fi

    # Fill in your application code here.

    echo "* Quantitate with extra flags: [${extraFlags}]..."
    set -x
    rsem-calculate-expression --bam --estimate-rspd --calc-ci --seed ${rnd_seed} -p ${nthreads} \
        --no-bam-output --ci-memory 30000 ${extraFlags} ${bam_fn}.bam ${index_prefix} ${bam_fn}_rsem
    set +x

    echo `ls ${bam_fn}*`

    echo "* Upload results..."
    rsem_gene_results=$(dx upload ${bam_fn}_rsem.genes.results   --property SW="$versions" --brief)
    rsem_iso_results=$(dx upload ${bam_fn}_rsem.isoforms.results --property SW="$versions" --brief)

    dx-jobutil-add-output rsem_gene_results "$rsem_gene_results" --class=file
    dx-jobutil-add-output rsem_iso_results "$rsem_iso_results" --class=file
    echo "* Finished."
}
