#!/bin/bash
# rampage-peaks.sh

main() {
    echo "* Installing python-dev cython python-scipy python-networkx..."
    sudo apt-get install -y gcc python-dev cython python-scipy python-networkx >> install.log 2>&1
    echo "* Installing pysam..."
    sudo easy_install pysam >> install.log 2>&1
    echo "* Installing grit..."
    wget https://github.com/nboley/grit/archive/2.0.5beta4.tar.gz -O grit.tgz
    mkdir grit_local
    tar -xzf grit.tgz -C grit_local --strip-components=1
    cd grit_local
    sudo python setup.py install >> ../install.log 2>&1
    cd ..

    # If available, will print tool versions to stderr and json string to stdout
    versions=''
    if [ -f /usr/bin/tool_versions.py ]; then 
        versions=`tool_versions.py --dxjson dnanexus-executable.json`
    fi

    echo "* Value of rampage bam: '$rampage_marked_bam'"
    echo "* Value of control bam: '$control_bam'"
    echo "* Value of annotation:  '$gene_annotation'"
    echo "* Value of chrom_sizes: '$chrom_sizes'"
    echo "* Value of nthreads:    '$nthreads'"

    echo "* Download files..."
    bam_root=`dx describe "$rampage_marked_bam" --name`
    bam_root=${bam_root%.bam}
    bam_root=${bam_root%_rampage_star_marked}
    dx download "$rampage_marked_bam" -o ${bam_root}.bam
    echo "* Bam file: '${bam_root}.bam'"

    control_root=`dx describe "$control_bam" --name`
    control_root=${control_root%.bam}
    dx download "$control_bam" -o "$control_root".bam
    echo "* Control bam file: '${control_root}.bam'"

    annotation_root=`dx describe "$gene_annotation" --name`
    annotation_root=${annotation_root%.gtf.gz}
    dx download "$gene_annotation" -o "$annotation_root".gtf.gz
    echo "* Annotation file: '${annotation_root}.gtf.gz'"
    
    dx download "$chrom_sizes" -o chrom.sizes

    # DX/ENCODE independent script is found in resources/usr/bin
    echo "* ===== Calling DNAnexus and ENCODE independent script... ====="
    set -x
    ram-peaks.sh ${bam_root}.bam ${control_root}.bam ${annotation_root}.gtf.gz chrom.sizes $nthreads 
    set +x
    echo "* ===== Returned from dnanexus and encodeD independent script ====="
    peaks_root=${bam_root}_rampage_peaks
    
    echo "* Upload results..."
    rampage_peaks_bed=$(dx upload ${peaks_root}.bed.gz --property SW="$versions" --brief)
    rampage_peaks_bb=$(dx upload ${peaks_root}.bb      --property SW="$versions" --brief)
    rampage_peaks_gff=$(dx upload ${peaks_root}.gff.gz --property SW="$versions" --brief)
    rampage_peak_quants=$(dx upload ${peaks_root}_quant.tsv --property SW="$versions" --brief)

    dx-jobutil-add-output rampage_peaks_bed "$rampage_peaks_bed" --class=file
    dx-jobutil-add-output rampage_peaks_bb "$rampage_peaks_bb" --class=file
    dx-jobutil-add-output rampage_peaks_gff "$rampage_peaks_gff" --class=file
    dx-jobutil-add-output rampage_peak_quants "$rampage_peak_quants" --class=file
    
    echo "* Finished."
}
