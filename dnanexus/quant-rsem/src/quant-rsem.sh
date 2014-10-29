#!/bin/bash
# quant-rsem 0.0.2

main() {
    # Now in resources/usr/bin
    # echo "* Dowload and install RSEM..."
    # git clone https://github.com/bli25wisc/RSEM.git
    # (cd RSEM; git checkout 92b24279a3ecc72946e7e7c23149ad0d181f373a)
    # #### get correct commit bit from submodule
    # (cd RSEM; make)

    echo "*****"
    echo "* Running: quant_rsem.sh"
    echo "* RSEM version: "`rsem-calculate-expression --version | awk '{print $5}'`
    echo "*****"

    echo "* Value of annotation_bam: '$annotation_bam'"
    echo "* Value of rsem_index: '$rsem_index'"
    echo "* Value of paired: '$paired'"
    #echo "* Value of stranded: '$stranded'"
    echo "* Random number seed: '$rnd_seed'"
    echo "* Value of nthreads: '$nthreads'"

    echo "* Download files..."
    bam_fn=`dx describe "$annotation_bam" --name`
    bam_fn=${bam_fn%.bam}
    dx download "$annotation_bam" -o ${bam_fn}.bam
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

    # TODO: should this be moved to align-star scripts?
    #if [ "$paired" == "true" ]; then
    #    echo "* Sorting for paired-end bam..."
    #    # paired-end data, merge mates into one line before sorting, and un-merge after sorting
    #    cat <( samtools view -H ${bam_fn}.bam ) \
    #        <( samtools view -@ ${nthreads} ${bam_fn}.bam | awk '{printf $0 " "; getline; print}' | \
    #            sort -S 60G -T ./ | tr ' ' '\n' ) | \
    #        samtools view -@ ${nthreads} -bS - > ${bam_fn}_sorted.bam
    #else
    #    echo "* Sorting for single-end bam..."
    #    cat <( samtools view -H ${bam_fn}.bam ) \
    #        <( samtools view -@ ${nthreads} ${bam_fn}.bam | sort -S 60G -T ./ ) | \
    #        samtools view -@ ${nthreads} -bS - > ${bam_fn}_sorted.bam
    #fi

    echo "* Quantitate with extra flags: [${extraFlags}]..."
    rsem-calculate-expression --bam --estimate-rspd --calc-ci --seed ${rnd_seed} -p ${nthreads} \
        --no-bam-output --ci-memory 30000 ${extraFlags} ${bam_fn}.bam ${index_prefix} ${bam_fn}_rsem

    echo `ls ${bam_fn}*`

    echo "* Upload results..."
    genomic_quant=$(dx upload ${bam_fn}_rsem.genes.results --brief)
    transcript_quant=$(dx upload ${bam_fn}_rsem.isoforms.results --brief)

    dx-jobutil-add-output genomic_quant "$genomic_quant" --class=file
    dx-jobutil-add-output transcript_quant "$transcript_quant" --class=file
    echo "* Finished."
}
