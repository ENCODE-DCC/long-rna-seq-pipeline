#!/bin/bash
# small-rna-signals.sh

main() {
    # Now in resources/usr/bin
    #echo "* Download and install STAR..."
    #git clone https://github.com/alexdobin/STAR
    #(cd STAR; git checkout tags/STAR_2.4.2a)
    #(cd STAR; make)
    #wget https://github.com/ENCODE-DCC/kentUtils/archive/v302.1.0.tar.gz

    # If available, will print tool versions to stderr and json string to stdout
    versions=''
    if [ -f /usr/bin/tool_versions.py ]; then 
        versions=`tool_versions.py --dxjson dnanexus-executable.json`
    fi

    echo "Value of srna_bam: '$srna_bam'"
    echo "Value of chrom_sizes: '$chrom_sizes'"

    echo "* Download files..."
    bam_root=`dx describe "$srna_bam" --name`
    bam_root=${bam_root%.bam}
    echo "* Bam file: '"$bam_root".bam'"
    dx download "$srna_bam" -o "$bam_root".bam

    dx download "$chrom_sizes" -o chrom.sizes

    # DX/ENCODE independent script is found in resources/usr/bin
    echo "* ===== Calling DNAnexus and ENCODE independent script... ====="
    set -x
    srna_signals.sh ${bam_root}.bam chrom.sizes
    set +x
    echo "* ===== Returned from dnanexus and encodeD independent script ====="

    echo "* Upload results..."
    all_minus_bw=$(dx upload ${bam_root}_minusAll.bw     --property SW="$versions" --brief)
    all_plus_bw=$(dx upload ${bam_root}_plusAll.bw       --property SW="$versions" --brief)
    unique_minus_bw=$(dx upload ${bam_root}_minusUniq.bw --property SW="$versions" --brief)
    unique_plus_bw=$(dx upload ${bam_root}_plusUniq.bw   --property SW="$versions" --brief)

    dx-jobutil-add-output all_minus_bw "$all_minus_bw" --class=file
    dx-jobutil-add-output all_plus_bw "$all_plus_bw" --class=file
    dx-jobutil-add-output unique_minus_bw "$unique_minus_bw" --class=file
    dx-jobutil-add-output unique_plus_bw "$unique_plus_bw" --class=file

    echo "* Finished."
}
