#!/bin/bash
# rampage-idr.sh

script_name="rampage-idr.sh"
script_ver="1.0.1"

main() {
    echo "* Installing Anaconda3 (python3.4.3, numpy-1.9.2 matplotlib-1.4.3..."
    set -x
    wget https://3230d63b5fc54e62148e-c95ac804525aac4b6dba79b00b39d1d3.ssl.cf1.rackcdn.com/Anaconda3-2.2.0-Linux-x86_64.sh >> ../install.log 2>&1
    bash Anaconda3-2.2.0-Linux-x86_64.sh -b
    ana_bin=`echo ~/anaconda3/bin`
    # python symlink will interfere with python2.7
    rm ${ana_bin}/python
    #ls ${ana_bin}
    export PATH=${ana_bin}:$PATH
    python3 -V 2>&1 | tee -a install.log
    #sudo python3 -V 2>&1 | tee -a install.log
    #sudo ${ana_bin}/python3 -V 2>&1 | tee -a install.log
    set +x

    echo "* Installing idr..."
    set -x
    wget https://github.com/nboley/idr/archive/2.0.0.tar.gz -O idr.tgz >> ../install.log 2>&1
    #wget https://github.com/nboley/idr/archive/2.0.0beta5.tar.gz -O idr.tgz >> ../install.log 2>&1
    mkdir idr
    tar -xzf idr.tgz -C idr --strip-components=1
    cd idr
    # sudo does not see python3 so it requires ana_bin path
    sudo ${ana_bin}/python3 setup.py install >> ../install.log 2>&1
    cd ..
    set +x
    
    # If available, will print tool versions to stderr and json string to stdout
    versions=''
    if [ -f /usr/bin/tool_versions.py ]; then 
        versions=`tool_versions.py --applet $script_name --appver $script_ver`
    fi

    echo "Value of peaks_a:     '$peaks_a'"
    echo "Value of peaks_b:     '$peaks_b'"
    echo "Value of chrom_sizes: '$chrom_sizes'"

    echo "* Download files..."
    peaks_a_fn=`dx describe "$peaks_a" --name`
    peaks_a_fn=${peaks_a_fn%_rampage_peaks.bed}
    peaks_a_fn=${peaks_a_fn%.bed}
    dx download "$peaks_a" -o peaks_a.bed
    echo "* First bed: '"$peaks_a_fn".bed'"

    peaks_b_fn=`dx describe "$peaks_b" --name`
    peaks_b_fn=${peaks_b_fn%_rampage_peaks.bed}
    peaks_b_fn=${peaks_b_fn%.bed}
    dx download "$peaks_b" -o peaks_b.bed
    echo "* Second bed: '"$peaks_b_fn".bed'"

    dx download "$chrom_sizes" -o chromSizes.txt

    idr_root=${peaks_a_fn}_${peaks_b_fn}_idr
    echo "* Rampage IDR root: '"$idr_root"'"

    echo "* Removing any spike-ins from bed files..."
    set -x
    grep "^chr" peaks_a.bed > peaks_a_clean.bed
    grep "^chr" peaks_b.bed > peaks_b_clean.bed
    set -x

    echo "* Running IDR..."
    set -x
    idr/bin/idr --input-file-type bed --rank 7 --plot --verbose --samples peaks_a_clean.bed peaks_b_clean.bed 2>&1 | tee idr_summary.txt
    sort -k1,1 -k2,2n < idrValues.txt > ${idr_root}.bed
    mv idrValues.txt.png ${idr_root}.png
    set +x

    echo "* Converting bed to bigBed..."
    set -x
    bedToBigBed ${idr_root}.bed -type=bed6+ -as=/usr/bin/rampage_idr_peaks.as chromSizes.txt ${idr_root}.bb
    set +x

    echo "* Prepare metadata..."
    meta=''
    if [ -f /usr/bin/qc_metrics.py ]; then
        meta=`qc_metrics.py -n IDR_summary -f idr_summary.txt`
    fi
    
    echo "* Upload results..."
    # NOTE: adding meta 'details' ensures json is valid.  But details are not updatable so rely on QC property
    details=`echo { $meta }`
    rampage_idr_bed=$(dx upload ${idr_root}.bed --details="$details" --property QC="$meta" --property SW="$versions" --brief)
    rampage_idr_bb=$(dx upload ${idr_root}.bb   --details="$details" --property QC="$meta" --property SW="$versions" --brief)
    rampage_idr_png=$(dx upload ${idr_root}.png --details="$details" --property QC="$meta" --property SW="$versions" --brief)

    dx-jobutil-add-output rampage_idr_bed "$rampage_idr_bed" --class=file
    dx-jobutil-add-output rampage_idr_bb "$rampage_idr_bb" --class=file
    dx-jobutil-add-output rampage_idr_png "$rampage_idr_png" --class=file
    dx-jobutil-add-output metadata "$meta" --class=string

    echo "* Finished."
}
