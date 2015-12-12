#!/bin/bash
# bam-to-bigwig-stranded.sh

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

    echo "Value of bam_file: '$bam_file'"
    echo "Value of chrom_sizes: '$chrom_sizes'"

    echo "* Download files..."
    bam_root=`dx describe "$bam_file" --name`
    bam_root=${bam_root%.bam}
    echo "* Bam file: '"$bam_root".bam'"
    dx download "$bam_file" -o "$bam_root".bam

    dx download "$chrom_sizes" -o chromSizes.txt

    # DX/ENCODE independent script is found in resources/usr/bin
    echo "* ===== Calling DNAnexus and ENCODE independent script... ====="
    set -x
    lrna-bam-to-stranded-signals.sh ${bam_root}.bam chromSizes.txt
    set +x
    echo "* ===== Returned from dnanexus and encodeD independent script ====="

    echo "* Upload results..."
    minus_all_bw=$(dx upload ${bam_root}_minusAll.bw   --property SW="$versions" --brief)
    minus_uniq_bw=$(dx upload ${bam_root}_minusUniq.bw --property SW="$versions" --brief)
    plus_all_bw=$(dx upload ${bam_root}_plusAll.bw     --property SW="$versions" --brief)
    plus_uniq_bw=$(dx upload ${bam_root}_plusUniq.bw   --property SW="$versions" --brief)

    dx-jobutil-add-output minus_all_bw "$minus_all_bw" --class=file
    dx-jobutil-add-output minus_uniq_bw "$minus_uniq_bw" --class=file
    dx-jobutil-add-output plus_all_bw "$plus_all_bw" --class=file
    dx-jobutil-add-output plus_uniq_bw "$plus_uniq_bw" --class=file

    echo "* Finished."
}
