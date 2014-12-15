#!/bin/bash
# rampage-peaks 0.0.1

main() {
    echo "* Install grit depencies..."
    #sudo apt-get install -y python-support python-numpy libamd2.2.0 libumfpack5.4.0 python-scipy python-networkx
    sudo apt-get install -y gcc python-dev cython python-scipy python-networkx
    echo "* Install pysam..."
    #pip install pysam
    sudo easy_install pysam
    echo "* Download and install grit..."
    wget https://github.com/nboley/grit/archive/2.0.0.tar.gz -O grit.tgz
    #mkdir grit
    #tar -xzf grit.tgz -C grit --strip-components=1
    #cd grit
    #sudo python setup.py install
    sudo easy_install grit.tgz
    echo "*****"
    echo `ls`
    ls
    echo "*****"
    #cd grit-2.0.0

    # Now in resources/usr/bin
    #wget https://github.com/ENCODE-DCC/kentUtils/archive/v302.1.0.tar.gz

    echo "*****"
    echo "* Running: rampage-peaks.sh [v0.0.1]"
    # TODO: call_peaks.py version?
    echo "* Running: grit:call_peaks.py [v2.0.0]"
    echo "* bedToBigBed version: "`bedToBigBed 2>&1 | grep "bedToBigBed v" | awk '{printf "v%s", $3}'`
    echo "* samtools version: "`samtools 2>&1 | grep Version | awk '{print $2}'`
    echo "*****"

    echo "Value of rampage bam: '$rampage_marked_bam'"
    echo "Value of control bam: '$control_bam'"
    echo "Value of annotation:  '$gene_annotation'"
    echo "Value of chrom_sizes: '$chrom_sizes'"
    echo "Value of as_file:     '$bed12_as_file'"
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
    dx download "$bed12_as_file" -o bed12.as

    peaks_root=${bam_fn}_rampage_peaks
    echo "* Rampage peaks root: '${peaks_root}'"
    
    ################ This step needs to be rewritten using call_peaks.py
    
    echo "* Indexing bams..."
    samtools index ${bam_fn}.bam 
    samtools index ${control_fn}.bam 

    echo "* Calling peaks..."
    python2.7 call_peaks.py --rampage-reads ${bam_fn}.bam --rnaseq-reads ${control_fn}.bam \
                                --reference ${annotation_fn}.gtf --exp-filter-fraction 0.05 \
                                --trim-fraction 0.01 --threads $nthreads --ucsc \
                                --outfname-type narrowPeak --outfname ${peaks_root}.bed \
                                --gene-regions-ofname ${peaks_root}_regions.bed
    echo `ls`

    #python2.7 bin/call_peaks.py --rampage-reads ${bam_fn}.bam --rnaseq-reads ${control_fn}.bam \
    #                            --reference ${annotation_fn}.gtf --exp-filter-fraction 0.05 \
    #                            --trim-fraction 0.01 --threads $nthreads --ucsc \
    #                            --outfname-type gff --outfname ${peaks_root}.gff
    touch ${peaks_root}.gff
 
    echo "* Converting narrowPeak bed to bigBed..."
    #bedToBigBed ${peaks_root}.bed -as=bed12.as chromSizes.txt ${peaks_root}.bb
    bedToBigBed ${peaks_root}.bed -type=bed6+3 chromSizes.txt ${peaks_root}.bb

    echo "* Upload results..."
    rampage_peaks_bed=$(dx upload ${peaks_root}.bed --brief)
    rampage_peaks_bb=$(dx upload ${peaks_root}.bb --brief)
    rampage_peaks_gff=$(dx upload ${peaks_root}.gff --brief)

    dx-jobutil-add-output rampage_peaks_bed "$rampage_peaks_bed" --class=file
    dx-jobutil-add-output rampage_peaks_bb "$rampage_peaks_bb" --class=file
    dx-jobutil-add-output rampage_peaks_gff "$rampage_peaks_gff" --class=file
    
    # temporary
    rampage_regions_bed=$(dx upload ${peaks_root}_regions.bed --brief)
    dx-jobutil-add-output rampage_regions_bed "$rampage_regions_bed" --class=file

    echo "* Finished."
}
