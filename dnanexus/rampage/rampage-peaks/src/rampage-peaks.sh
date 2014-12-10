#!/bin/bash
# rampage-peaks 0.0.1

main() {
    # Now in resources/usr/bin
    #echo "* Download and install STAR..."
    #git clone https://github.com/nboley/grit.git
    #cd grit
    #sudo apt-get install python-support python-numpy libamd2.2.0 libumfpack5.4.0 python-scipy python-networkx
    
    #(cd STAR; git checkout tags/STAR_2.4.0d)
    #(cd STAR; make)
    #wget https://github.com/ENCODE-DCC/kentUtils/archive/v302.1.0.tar.gz

    echo "*****"
    echo "* Running: rampage-peaks.sh [v0.0.1]"
    echo "* STAR version:     ["`STAR --version | awk '{print $1}' | cut -d _ -f 2-`"]"
    echo "* bedToBigBed version: "`bedToBigBed 2>&1 | grep "bedToBigBed v" | awk '{printf "v%s", $3}'`
    # TODO: rampagePeakCaller.py version?
    echo "* samtools version: "`samtools 2>&1 | grep Version | awk '{print $2}'`
    echo "*****"

    echo "Value of marked_bam:  '$marked_bam'"
    echo "Value of marked_bam:  '$control_bam'"
    echo "Value of annotation:  '$annotation_gtf'"
    echo "Value of chrom_sizes: '$chrom_sizes'"
    echo "* Value of nthreads:  '$nthreads'"

    echo "* Download files..."
    bam_fn=`dx describe "$marked_bam" --name`
    bam_fn=${bam_fn%.bam}
    dx download "$marked_bam" -o "$bam_fn".bam
    echo "* Bam file: '"$bam_fn".bam'"

    control_fn=`dx describe "$control_bam" --name`
    control_fn=${control_fn%.bam}
    dx download "$control_bam" -o "$control_fn".bam
    echo "* Control bam file: '"$control_fn".bam'"

    annotation_fn=`dx describe "$annotation_gtf" --name`
    annotation_fn=${annotation_fn%.gtf.gz}
    dx download "$annotation_gtf" -o "$annotation_fn".gtf.gz
    gunzip "$annotation_fn".gtf.gz
    echo "* Annotation file: '"$annotation_fn".gtf'"
    
    dx download "$chrom_sizes" -o chromSizes.txt

    ################ This step needs to be rewritten using call_peaks.py
    
    echo "* Inndexing bams..."
    samtools index ${bam_fn}.bam 
    samtools index ${control_fn}.bam 

    echo "* Calling peaks..."
    python2.7 /usr/bin/call_peaks.py --rampage-reads ${bam_fn}.bam --threads ${nthreads} \
                                     --rnaseq-reads ${control_fn}.bam --quiet \
                                     --reference ${annotation_fn}.gtf --exp-filter-fraction 0.05 \
                                     --trim-fraction 0.01 --out-fname ${bam_fn}_rampage_peaks.bed 

    echo "* Converting narrowPeak bed to bigBed..."
    bedToBigBed ${bam_fn}_rampage_peaks.bed -type=bed12 chromSizes.txt ${bam_fn}_rampage_peaks.bb

    echo "* Upload results..."
    rampage_peaks=$(dx upload ${bam_fn}_rampage_peaks.bed --brief)
    rampage_peaks_bb=$(dx upload ${bam_fn}_rampage_peaks.bb --brief)

    dx-jobutil-add-output rampage_peaks "$rampage_peaks" --class=file
    dx-jobutil-add-output rampage_peaks_bb "$rampage_peaks_bb" --class=file

    echo "* Finished."
}
