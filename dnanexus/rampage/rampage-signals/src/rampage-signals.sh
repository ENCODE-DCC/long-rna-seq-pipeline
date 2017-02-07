#!/bin/bash
# rampage-signals.sh

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
        versions=`tool_versions.py --dxjson dnanexus-executable.json`
    fi

    echo "Value of bam_file:    '$rampage_marked_bam'"
    echo "Value of chrom_sizes: '$chrom_sizes'"
    echo "Value of stranded:    '$stranded'"

    echo "* Download files..."
    bam_root=`dx describe "$rampage_marked_bam" --name`
    bam_root=${bam_root%.bam}
    bam_root=${bam_root%_star_marked}
    assay_type=${bam_root##*_}
    bam_root=${bam_root%_rampage}
    bam_root=${bam_root%_cage}
    dx download "$rampage_marked_bam" -o "$bam_root".bam
    echo "* Bam file: '${bam_root}.bam'"

    dx download "$chrom_sizes" -o chrom.sizes
    
    # DX/ENCODE independent script is found in resources/usr/bin
    signal_root=${bam_root}_${assay_type}_5p
    echo "* ===== Calling DNAnexus and ENCODE independent script... ====="
    set -x
    rampage_signal.sh ${bam_root}.bam chrom.sizes $signal_root $stranded 
    set +x
    echo "* ===== Returned from dnanexus and encodeD independent script ====="

    echo "* Upload results..."
    if [ "$stranded" == "true" ]; then
        all_minus_bw=$(dx upload ${signal_root}_minusAll.bw     --property SW="$versions" --brief)
        all_plus_bw=$(dx upload ${signal_root}_plusAll.bw       --property SW="$versions" --brief)
        unique_minus_bw=$(dx upload ${signal_root}_minusUniq.bw --property SW="$versions" --brief)
        unique_plus_bw=$(dx upload ${signal_root}_plusUniq.bw   --property SW="$versions" --brief)

        dx-jobutil-add-output all_minus_bw "$all_minus_bw" --class=file
        dx-jobutil-add-output all_plus_bw "$all_plus_bw" --class=file
        dx-jobutil-add-output unique_minus_bw "$unique_minus_bw" --class=file
        dx-jobutil-add-output unique_plus_bw "$unique_plus_bw" --class=file
    else
        all_bw=$(dx upload ${signal_root}_all.bw   --property SW="$versions" --brief)
        uniq_bw=$(dx upload ${signal_root}_uniq.bw --property SW="$versions" --brief)

        dx-jobutil-add-output all_bw "$all_bw" --class=file
        dx-jobutil-add-output uniq_bw "$uniq_bw" --class=file
    fi

    echo "* Finished."
}
