#!/bin/bash
# quant-rsem 1.0.1

main() {
    # Now in resources/usr/bin
    # echo "* Dowload and install RSEM..."
    # git clone https://github.com/bli25wisc/RSEM.git
    # (cd RSEM; git checkout tags/v1.2.19)
    # #### get correct commit bit from submodule
    # (cd RSEM; make)

    echo "*****"
    echo "* Running: quant_rsem.sh [v1.0.1]"
    echo "* RSEM version: "`rsem-calculate-expression --version | awk '{print $5}'`
    echo "*****"

    echo "* Value of annotation_bam: '$star_anno_bam'"
    echo "* Value of rsem_index: '$rsem_index'"
    echo "* Value of paired: '$paired'"
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
    if [ "$paired" == "true" ]; then
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
    rsem-calculate-expression --bam --estimate-rspd --calc-ci --seed ${rnd_seed} -p ${nthreads} \
        --no-bam-output --ci-memory 30000 ${extraFlags} ${bam_fn}.bam ${index_prefix} ${bam_fn}_rsem

    echo `ls ${bam_fn}*`

    echo "* Upload results..."
    rsem_gene_results=$(dx upload ${bam_fn}_rsem.genes.results --brief)
    rsem_iso_results=$(dx upload ${bam_fn}_rsem.isoforms.results --brief)

    dx-jobutil-add-output rsem_gene_results "$rsem_gene_results" --class=file
    dx-jobutil-add-output rsem_iso_results "$rsem_iso_results" --class=file
    echo "* Finished."
}
