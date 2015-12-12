#!/bin/bash
# bam-to-bigwig-unstranded.sh

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

    echo "* Value of bam_file: '$bam_file'"
    echo "* Value of chrom_sizes: '$chrom_sizes'"

    echo "* Download files..."
    bam_root=`dx describe "$bam_file" --name`
    bam_root=${bam_root%.bam}
    echo "* Bam file: '"$bam_root".bam'"
    dx download "$bam_file" -o "$bam_root".bam

    dx download "$chrom_sizes" -o chromSizes.txt

    # DX/ENCODE independent script is found in resources/usr/bin
    echo "* ===== Calling DNAnexus and ENCODE independent script... ====="
    set -x
    lrna-bam-to-unstranded-signals.sh ${bam_root}.bam chromSizes.txt
    set +x
    echo "* ===== Returned from dnanexus and encodeD independent script ====="

    echo "* Upload results..."
    all_bw=$(dx upload ${bam_root}_all.bw   --property SW="$versions" --brief)
    uniq_bw=$(dx upload ${bam_root}_uniq.bw --property SW="$versions" --brief)

    dx-jobutil-add-output all_bw "$all_bw" --class=file
    dx-jobutil-add-output uniq_bw "$uniq_bw" --class=file

    echo "* Finished."
}
