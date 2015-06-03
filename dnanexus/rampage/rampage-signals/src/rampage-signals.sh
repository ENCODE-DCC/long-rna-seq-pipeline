#!/bin/bash
# rampage-signals.sh

script_name="rampage-signals.sh"
script_ver="1.0.1"

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

    echo "Value of bam_file:    '$rampage_marked_bam'"
    echo "Value of chrom_sizes: '$chrom_sizes'"

    echo "* Download files..."
    bam_fn=`dx describe "$rampage_marked_bam" --name`
    bam_fn=${bam_fn%_rampage_star_marked.bam}
    bam_fn=${bam_fn%.bam}
    dx download "$rampage_marked_bam" -o "$bam_fn".bam
    echo "* Bam file: '${bam_fn}.bam'"

    dx download "$chrom_sizes" -o chromSizes.txt
    
    signal_root=${bam_fn}_rampage_5p
    echo "* Signal files root: '${signal_root}'"

    echo "* Make signals..."
    set -x
    mkdir -p Signal
    STAR --runMode inputAlignmentsFromBAM --inputBAMfile ${bam_fn}.bam --outWigType bedGraph read1_5p \
         --outWigStrand Stranded --outFileNamePrefix read1_5p. --outWigReferencesPrefix chr
    set +x

    echo "* Convert bedGraph to bigWigs..."
    set -x
    bedGraphToBigWig read1_5p.Signal.UniqueMultiple.str2.out.bg chromSizes.txt ${signal_root}_minusAll.bw
    bedGraphToBigWig read1_5p.Signal.Unique.str2.out.bg         chromSizes.txt ${signal_root}_minusUniq.bw
    bedGraphToBigWig read1_5p.Signal.UniqueMultiple.str1.out.bg chromSizes.txt ${signal_root}_plusAll.bw
    bedGraphToBigWig read1_5p.Signal.Unique.str1.out.bg         chromSizes.txt ${signal_root}_plusUniq.bw
    set +x
    echo `ls`

    echo "* Upload results..."
    all_minus_bw=$(dx upload ${signal_root}_minusAll.bw     --property SW="$versions" --brief)
    all_plus_bw=$(dx upload ${signal_root}_plusAll.bw       --property SW="$versions" --brief)
    unique_minus_bw=$(dx upload ${signal_root}_minusUniq.bw --property SW="$versions" --brief)
    unique_plus_bw=$(dx upload ${signal_root}_plusUniq.bw   --property SW="$versions" --brief)

    dx-jobutil-add-output all_minus_bw "$all_minus_bw" --class=file
    dx-jobutil-add-output all_plus_bw "$all_plus_bw" --class=file
    dx-jobutil-add-output unique_minus_bw "$unique_minus_bw" --class=file
    dx-jobutil-add-output unique_plus_bw "$unique_plus_bw" --class=file

    echo "* Finished."
}
