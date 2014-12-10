#!/bin/bash
# rampage-signals 0.0.1

main() {
    # Now in resources/usr/bin
    #echo "* Download and install STAR..."
    #git clone https://github.com/alexdobin/STAR
    #(cd STAR; git checkout tags/STAR_2.4.0g1)
    #(cd STAR; make)
    #wget https://github.com/ENCODE-DCC/kentUtils/archive/v302.1.0.tar.gz

    echo "*****"
    echo "* Running: rampage-signals.sh [v0.0.1]"
    echo "* STAR version:     ["`STAR --version | awk '{print $1}' | cut -d _ -f 2-`"]"
    echo "* bedGraphToBigWig version: "`bedGraphToBigWig 2>&1 | grep "bedGraphToBigWig v" | awk '{print $2$3}'`
    echo "*****"

    echo "Value of bam_file: '$bam_file'"
    echo "Value of chrom_sizes: '$chrom_sizes'"

    echo "* Download files..."
    bam_fn=`dx describe "$bam_file" --name`
    bam_fn=${bam_fn%.bam}
    echo "* Bam file: '"$bam_fn".bam'"
    dx download "$bam_file" -o "$bam_fn".bam

    dx download "$chrom_sizes" -o chromSizes.txt

    echo "* Make signals..."
    mkdir -p Signal
    STAR --runMode inputAlignmentsFromBAM --inputBAMfile ${bam_fn}.bam --outWigType bedGraph read1_5p \
         --outWigStrand Stranded --outFileNamePrefix read1_5p. --outWigReferencesPrefix chr

    echo "* Convert bedGraph to bigWigs..."
    bedGraphToBigWig read1_5p.Signal.UniqueMultiple.str2.out.bg chromSizes.txt ${bam_fn}_5p_minusAll.bw
    bedGraphToBigWig read1_5p.Signal.Unique.str2.out.bg         chromSizes.txt ${bam_fn}_5p_minusUniq.bw
    bedGraphToBigWig read1_5p.Signal.UniqueMultiple.str1.out.bg chromSizes.txt ${bam_fn}_5p_plusAll.bw
    bedGraphToBigWig read1_5p.Signal.Unique.str1.out.bg         chromSizes.txt ${bam_fn}_5p_plusUniq.bw
    echo `ls`

    echo "* Upload results..."
    all_minus_bw=$(dx upload ${bam_fn}_5p_minusAll.bw --brief)
    all_plus_bw=$(dx upload ${bam_fn}_5p_plusAll.bw --brief)
    unique_minus_bw=$(dx upload ${bam_fn}_5p_minusUniq.bw --brief)
    unique_plus_bw=$(dx upload ${bam_fn}_5p_plusUniq.bw --brief)

    dx-jobutil-add-output all_minus_bw "$all_minus_bw" --class=file
    dx-jobutil-add-output all_plus_bw "$all_plus_bw" --class=file
    dx-jobutil-add-output unique_minus_bw "$unique_minus_bw" --class=file
    dx-jobutil-add-output unique_plus_bw "$unique_plus_bw" --class=file

    #echo "* Temporary uploads..."
    # temprary for comparison only!
    mv read1_5p.Signal.UniqueMultiple.str2.out.bg ${bam_fn}_5p_minusAll.bg
    mv read1_5p.Signal.Unique.str2.out.bg         ${bam_fn}_5p_minusUniq.bg
    mv read1_5p.Signal.UniqueMultiple.str1.out.bg ${bam_fn}_5p_plusAll.bg
    mv read1_5p.Signal.Unique.str1.out.bg         ${bam_fn}_5p_plusUniq.bg
    all_minus_bg=$(dx upload ${bam_fn}_5p_minusAll.bg --brief)
    all_plus_bg=$(dx upload ${bam_fn}_5p_plusAll.bg --brief)
    unique_minus_bg=$(dx upload ${bam_fn}_5p_minusUniq.bg --brief)
    unique_plus_bg=$(dx upload ${bam_fn}_5p_plusUniq.bg --brief)
    dx-jobutil-add-output all_minus_bg "$all_minus_bg" --class=file
    dx-jobutil-add-output all_plus_bg "$all_plus_bg" --class=file
    dx-jobutil-add-output unique_minus_bg "$unique_minus_bg" --class=file
    dx-jobutil-add-output unique_plus_bg "$unique_plus_bg" --class=file
    echo "* Finished."
}
