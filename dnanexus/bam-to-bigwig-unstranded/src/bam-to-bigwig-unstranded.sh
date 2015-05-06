#!/bin/bash
# bam-to-bigwig-unstranded.sh

script_name="bam-to-bigwig-unstranded.sh"
script_ver="2.0.1"

main() {
    # Now in resources/usr/bin
    #echo "* Download and install STAR..."
    #git clone https://github.com/alexdobin/STAR
    #(cd STAR; git checkout tags/STAR_2.4.0k)
    #(cd STAR; make)
    #wget https://github.com/ENCODE-DCC/kentUtils/archive/v302.1.0.tar.gz

    # If available, will print tool versions to stderr and json string to stdout
    versions=''
    if [ -f /usr/bin/tool_versions.py ]; then 
        versions=`tool_versions.py --applet $script_name --appver $script_ver`
    fi
    #echo "*****"
    #echo "* Running: bam-to-bigwig-unstranded.sh [v2.0.0]"
    #echo "* STAR version:     ["`STAR --version | awk '{print $1}' | cut -d _ -f 2-`"]"
    #echo "* bedGraphToBigWig version: "`bedGraphToBigWig 2>&1 | grep "bedGraphToBigWig v" | awk '{print $2$3}'`
    #echo "*****"

    echo "* Value of bam_file: '$bam_file'"
    echo "* Value of chrom_sizes: '$chrom_sizes'"

    # The following line(s) use the dx command-line tool to download your file
    # inputs to the local file system using variable names for the filenames. To
    # recover the original filenames, you can use the output of "dx describe
    # "$variable" --name".

    echo "* Download files..."
    bam_fn=`dx describe "$bam_file" --name`
    bam_fn=${bam_fn%.bam}
    echo "* Bam file: '"$bam_fn".bam'"
    dx download "$bam_file" -o "$bam_fn".bam

    dx download "$chrom_sizes" -o chromSizes.txt

    echo "* Make signals..."
    set -x
    mkdir -p Signal
    STAR --runMode inputAlignmentsFromBAM --inputBAMfile ${bam_fn}.bam --outWigType bedGraph \
         --outWigStrand Unstranded --outFileNamePrefix ./Signal/ --outWigReferencesPrefix chr
    mv Signal/Signal*bg .
    set +x
    echo `ls -l`

    echo "* Convert bedGraph to bigWigs..."
    set -x
    bedGraphToBigWig Signal.UniqueMultiple.str1.out.bg chromSizes.txt ${bam_fn}_all.bw
    bedGraphToBigWig Signal.Unique.str1.out.bg         chromSizes.txt ${bam_fn}_uniq.bw
    set +x

    echo "* Upload results..."
    all_bw=$(dx upload ${bam_fn}_all.bw   --property SW="$versions" --brief)
    uniq_bw=$(dx upload ${bam_fn}_uniq.bw --property SW="$versions" --brief)

    dx-jobutil-add-output all_bw "$all_bw" --class=file
    dx-jobutil-add-output uniq_bw "$uniq_bw" --class=file

    echo "* Finished."
}
