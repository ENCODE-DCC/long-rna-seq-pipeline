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
    echo "* Value of rsem_index:     '$rsem_index'"
    echo "* Value of paired:         '$paired_end'"
    echo "* Random number seed:      '$rnd_seed'"
    echo "* Value of nthreads:       '$nthreads'"

    echo "* Download files..."
    bam_root=`dx describe "$star_anno_bam" --name`
    bam_root=${bam_root%_star_anno.bam}
    bam_root=${bam_root%.bam}
    dx download "$star_anno_bam" -o ${bam_root}.bam
    dx download "$rsem_index" -o rsem_index.tgz

    # DX/ENCODE independent script is found in resources/usr/bin
    echo "* ===== Calling DNAnexus and ENCODE independent script... ====="
    set -x
    lrna-rsem-quantification.sh rsem_index.tgz ${bam_root}.bam $paired_end $rnd_seed $nthreads
    set +x
    echo "* ===== Returned from dnanexus and encodeD independent script ====="

    echo "* Upload results..."
    rsem_gene_results=$(dx upload ${bam_root}_rsem.genes.results   --property SW="$versions" --brief)
    rsem_iso_results=$(dx upload ${bam_root}_rsem.isoforms.results --property SW="$versions" --brief)

    dx-jobutil-add-output rsem_gene_results "$rsem_gene_results" --class=file
    dx-jobutil-add-output rsem_iso_results "$rsem_iso_results" --class=file
    echo "* Finished."
}
