#!/bin/bash
# rampage-peaks 1.0.0

main() {
    echo "* Installing python-dev cython python-scipy python-networkx..."
    #sudo apt-get install -y python-support python-numpy libamd2.2.0 libumfpack5.4.0 python-scipy python-networkx
    sudo apt-get install -y gcc python-dev cython python-scipy python-networkx >> install.log 2>&1
    echo "* Installing pysam..."
    #pip install pysam
    sudo easy_install pysam >> install.log 2>&1
    echo "* Installing grit..."
    wget https://github.com/nboley/grit/archive/2.0.2.tar.gz -O grit.tgz
    mkdir grit_local
    tar -xzf grit.tgz -C grit_local --strip-components=1
    #git clone https://github.com/nboley/grit.git
    #mv grit grit_local
    cd grit_local
    sudo python setup.py install >> ../install.log 2>&1
    #echo "***** ls ."
    #ls
    #echo "***** ls /usr/local/bin"
    #ls /usr/local/bin
    cd ..

    # Now in resources/usr/bin
    #wget https://github.com/ENCODE-DCC/kentUtils/archive/v302.1.0.tar.gz

    echo "*****"
    echo "* Running: rampage-peaks.sh [v1.0.0]"
    echo "* Running: grit:call_peaks version: "`call_peaks --version 2>&1 | grep call_peaks | awk '{print $3}'`
    echo "* bedToBigBed version: "`bedToBigBed 2>&1 | grep "bedToBigBed v" | awk '{printf "v%s", $3}'`
    echo "* samtools version: "`samtools 2>&1 | grep Version | awk '{print $2}'`
    echo "*****"

    echo "Value of rampage bam: '$rampage_marked_bam'"
    echo "Value of control bam: '$control_bam'"
    echo "Value of annotation:  '$gene_annotation'"
    echo "Value of chrom_sizes: '$chrom_sizes'"
    echo "* Value of nthreads:  '$nthreads'"

    echo "* Download files..."
    bam_fn=`dx describe "$rampage_marked_bam" --name`
    bam_fn=${bam_fn%_rampage_star_marked.bam}
    bam_fn=${bam_fn%.bam}
    dx download "$rampage_marked_bam" -o "$bam_fn".bam
    echo "* Bam file: '${bam_fn}.bam'"

    control_fn=`dx describe "$control_bam" --name`
    control_fn=${control_fn%.bam}
    dx download "$control_bam" -o "$control_fn".bam
    echo "* Control bam file: '${control_fn}.bam'"

    annotation_fn=`dx describe "$gene_annotation" --name`
    annotation_fn=${annotation_fn%.gtf.gz}
    dx download "$gene_annotation" -o "$annotation_fn".gtf.gz
    gunzip "$annotation_fn".gtf.gz
    echo "* Annotation file: '${annotation_fn}.gtf'"
    
    dx download "$chrom_sizes" -o chromSizes.txt

    peaks_root=${bam_fn}_rampage_peaks
    echo "* Rampage peaks root: '${peaks_root}'"
    
    # TODO:
    # Should we remove the ERCC spike-in peaks first?
    echo "* Indexing bams..."
    samtools index ${bam_fn}.bam 
    samtools index ${control_fn}.bam 

    echo "* Calling peaks..."
    call_peaks --rampage-reads ${bam_fn}.bam --rnaseq-reads ${control_fn}.bam --threads $nthreads \
               --reference ${annotation_fn}.gtf --exp-filter-fraction 0.05 --trim-fraction 0.01 \
                --ucsc --outfname ${peaks_root}.gff --outfname-type gff \
               --bed-peaks-ofname ${peaks_root}.bed
    echo `ls`
 
    echo "* Converting bed to bigBed..."
    grep "^chr" ${peaks_root}.bed | sort -k1,1 -k2,2n > peaks_polished.bed
    #cp peaks_polished.bed ${peaks_root}.bb
    bedToBigBed peaks_polished.bed -type=bed6+ -as=/usr/bin/tss_peak.as chromSizes.txt ${peaks_root}.bb

    echo "* Upload results..."
    rampage_peaks_bed=$(dx upload ${peaks_root}.bed --brief)
    rampage_peaks_bb=$(dx upload ${peaks_root}.bb --brief)
    rampage_peaks_gff=$(dx upload ${peaks_root}.gff --brief)

    dx-jobutil-add-output rampage_peaks_bed "$rampage_peaks_bed" --class=file
    dx-jobutil-add-output rampage_peaks_bb "$rampage_peaks_bb" --class=file
    dx-jobutil-add-output rampage_peaks_gff "$rampage_peaks_gff" --class=file
    
    echo "* Finished."
}
