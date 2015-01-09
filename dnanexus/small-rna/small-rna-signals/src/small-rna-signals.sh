#!/bin/bash
# small-rna-signals 1.0.0

main() {
    # Now in resources/usr/bin
    #echo "* Download and install STAR..."
    #git clone https://github.com/alexdobin/STAR
    #(cd STAR; git checkout tags/STAR_2.4.0g1)
    #(cd STAR; make)
    #wget https://github.com/ENCODE-DCC/kentUtils/archive/v302.1.0.tar.gz

    echo "*****"
    echo "* Running: small-rna-signals.sh [v1.0.0]"
    echo "* STAR version:     ["`STAR --version | awk '{print $1}' | cut -d _ -f 2-`"]"
    echo "* bedGraphToBigWig version: "`bedGraphToBigWig 2>&1 | grep "bedGraphToBigWig v" | awk '{print $2$3}'`
    echo "*****"

    echo "Value of genome_bam: '$genome_bam'"
    echo "Value of chrom_sizes: '$chrom_sizes'"

    echo "* Download files..."
    bam_fn=`dx describe "$genome_bam" --name`
    bam_fn=${bam_fn%.bam}
    echo "* Bam file: '"$bam_fn".bam'"
    dx download "$genome_bam" -o "$bam_fn".bam

    dx download "$chrom_sizes" -o chromSizes.txt

    echo "* Make signals..."
    mkdir -p Signal
    STAR --runMode inputAlignmentsFromBAM --inputBAMfile ${bam_fn}.bam --outWigType bedGraph \
         --outWigStrand Stranded --outWigReferencesPrefix chr

    echo "* Convert bedGraph to bigWigs..."
    bedGraphToBigWig Signal.UniqueMultiple.str2.out.bg chromSizes.txt ${bam_fn}_small_minusAll.bw
    bedGraphToBigWig Signal.Unique.str2.out.bg         chromSizes.txt ${bam_fn}_small_minusUniq.bw
    bedGraphToBigWig Signal.UniqueMultiple.str1.out.bg chromSizes.txt ${bam_fn}_small_plusAll.bw
    bedGraphToBigWig Signal.Unique.str1.out.bg         chromSizes.txt ${bam_fn}_small_plusUniq.bw
    echo `ls`

    echo "* Upload results..."
    all_minus_bw=$(dx upload ${bam_fn}_small_minusAll.bw --brief)
    all_plus_bw=$(dx upload ${bam_fn}_small_plusAll.bw --brief)
    unique_minus_bw=$(dx upload ${bam_fn}_small_minusUniq.bw --brief)
    unique_plus_bw=$(dx upload ${bam_fn}_small_plusUniq.bw --brief)

    dx-jobutil-add-output all_minus_bw "$all_minus_bw" --class=file
    dx-jobutil-add-output all_plus_bw "$all_plus_bw" --class=file
    dx-jobutil-add-output unique_minus_bw "$unique_minus_bw" --class=file
    dx-jobutil-add-output unique_plus_bw "$unique_plus_bw" --class=file

    #echo "* Temporary uploads..."
    # temprary for comparison only!
    mv Signal.UniqueMultiple.str2.out.bg ${bam_fn}_small_minusAll.bg
    mv Signal.Unique.str2.out.bg         ${bam_fn}_small_minusUniq.bg
    mv Signal.UniqueMultiple.str1.out.bg ${bam_fn}_small_plusAll.bg
    mv Signal.Unique.str1.out.bg         ${bam_fn}_small_plusUniq.bg
    all_minus_bg=$(dx upload ${bam_fn}_small_minusAll.bg --brief)
    all_plus_bg=$(dx upload ${bam_fn}_small_plusAll.bg --brief)
    unique_minus_bg=$(dx upload ${bam_fn}_small_minusUniq.bg --brief)
    unique_plus_bg=$(dx upload ${bam_fn}_small_plusUniq.bg --brief)
    dx-jobutil-add-output all_minus_bg "$all_minus_bg" --class=file
    dx-jobutil-add-output all_plus_bg "$all_plus_bg" --class=file
    dx-jobutil-add-output unique_minus_bg "$unique_minus_bg" --class=file
    dx-jobutil-add-output unique_plus_bg "$unique_plus_bg" --class=file
    echo "* Finished."
}
