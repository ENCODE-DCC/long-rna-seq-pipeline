#!/bin/bash
# bam-to-bigwig 1.0.0

main() {
    # Now in resources/usr/bin
    #echo "* Download and install STAR..."
    #git clone https://github.com/alexdobin/STAR
    #(cd STAR; git checkout tags/STAR_2.4.0d)
    #(cd STAR; make)
    #wget https://github.com/ENCODE-DCC/kentUtils/archive/v302.1.0.tar.gz

    echo "*****"
    echo "* Running: bam-to-bigwig-stranded.sh [v1.0.0]"
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
    STAR --runMode inputAlignmentsFromBAM --inputBAMfile ${bam_fn}.bam --outWigType bedGraph \
         --outWigStrand Stranded --outFileNamePrefix ./Signal/ --outWigReferencesPrefix chr
    mv Signal/Signal*bg .

    echo "* Convert bedGraph to bigWigs..."
    bedGraphToBigWig Signal.UniqueMultiple.str1.out.bg chromSizes.txt ${bam_fn}_minusAll.bw
    bedGraphToBigWig Signal.Unique.str1.out.bg         chromSizes.txt ${bam_fn}_minusUniq.bw
    bedGraphToBigWig Signal.UniqueMultiple.str2.out.bg chromSizes.txt ${bam_fn}_plusAll.bw
    bedGraphToBigWig Signal.Unique.str2.out.bg         chromSizes.txt ${bam_fn}_plusUniq.bw
    echo `ls`

    echo "* Upload results..."
    all_minus_bw=$(dx upload ${bam_fn}_minusAll.bw --brief)
    all_plus_bw=$(dx upload ${bam_fn}_plusAll.bw --brief)
    unique_minus_bw=$(dx upload ${bam_fn}_minusUniq.bw --brief)
    unique_plus_bw=$(dx upload ${bam_fn}_plusUniq.bw --brief)

    dx-jobutil-add-output all_minus_bw "$all_minus_bw" --class=file
    dx-jobutil-add-output all_plus_bw "$all_plus_bw" --class=file
    dx-jobutil-add-output unique_minus_bw "$unique_minus_bw" --class=file
    dx-jobutil-add-output unique_plus_bw "$unique_plus_bw" --class=file

    #echo "* Temporary uploads..."
    # temprary for comparison only!
    #mv Signal.UniqueMultiple.str1.out.bg ${bam_fn}_minusAll.bg
    #mv Signal.Unique.str1.out.bg         ${bam_fn}_minusUniq.bg
    #mv Signal.UniqueMultiple.str2.out.bg ${bam_fn}_plusAll.bg
    #mv Signal.Unique.str2.out.bg         ${bam_fn}_plusUniq.bg
    #minusAll_bg=$(dx upload ${bam_fn}_minusAll.bg --brief)
    #plusAll_bg=$(dx upload ${bam_fn}_plusAll.bg --brief)
    #minusUniq_bg=$(dx upload ${bam_fn}_minusUniq.bg --brief)
    #plusUniq_bg=$(dx upload ${bam_fn}_plusUniq.bg --brief)
    #dx-jobutil-add-output minusAll_bg "$minusAll_bg" --class=file
    #dx-jobutil-add-output plusAll_bg "$plusAll_bg" --class=file
    #dx-jobutil-add-output minusUniq_bg "$minusUniq_bg" --class=file
    #dx-jobutil-add-output plusUniq_bg "$plusUniq_bg" --class=file
    echo "* Finished."
}
